//! A per-kind entity table with soft removal.

use std::ops::{Index, IndexMut};

use super::id::Id;
use super::revision::Revision;

/// Dense storage for one entity kind. Removal is soft (the slot stays, flagged), so ids of
/// other entries stay valid; [`Table::compact`] reclaims removed slots and renumbers. A slot
/// released by `remove_releasing` (which nothing live refers to) is reused by a later
/// [`Table::push`].
#[derive(Debug, Clone)]
pub struct Table<T> {
    items: Vec<T>,
    stable_ids: Vec<u32>,
    revisions: Vec<Revision>,
    /// Per entry, how often what it is (not its meshes) was edited: what values computed
    /// from it go by (see [`Table::get_mut_keeping_content`]).
    contents: Vec<u32>,
    removed: Vec<bool>,
    active: usize,
    /// Released slots, for reuse (one revived since is skipped).
    free: Vec<u32>,
    released: Vec<bool>,
}

impl<T> Default for Table<T> {
    fn default() -> Self {
        Self {
            items: Vec::new(),
            stable_ids: Vec::new(),
            revisions: Vec::new(),
            contents: Vec::new(),
            removed: Vec::new(),
            active: 0,
            free: Vec::new(),
            released: Vec::new(),
        }
    }
}

impl<T> Table<T> {
    /// Adds `item` with persistent counter `stable_id`, in a released slot if there is one
    /// (its counters carry on, so values derived from what was there do not match).
    pub fn push(&mut self, item: T, stable_id: u32) -> Id<T> {
        while let Some(i) = self.free.pop() {
            let i = i as usize;
            if !std::mem::take(&mut self.released[i]) {
                continue;
            }
            self.items[i] = item;
            self.stable_ids[i] = stable_id;
            let revision = &mut self.revisions[i];
            (revision.own, revision.subtree) = (revision.own.wrapping_add(1), revision.subtree.wrapping_add(1));
            self.contents[i] = self.contents[i].wrapping_add(1);
            self.removed[i] = false;
            self.active += 1;
            return Id::from_index(i as u32);
        }
        let id = Id::from_index(u32::try_from(self.items.len()).expect("table overflow"));
        self.items.push(item);
        self.stable_ids.push(stable_id);
        self.revisions.push(Revision::default());
        self.contents.push(0);
        self.removed.push(false);
        self.released.push(false);
        self.active += 1;
        id
    }

    /// The entry, or `None` if `id` is out of range or removed.
    pub fn get(&self, id: Id<T>) -> Option<&T> {
        let i = id.index() as usize;
        (i < self.items.len() && !self.removed[i]).then(|| &self.items[i])
    }

    /// The entry for changing it (which counts as an edit: its own revision moves on), or
    /// `None` if `id` is out of range or removed.
    pub fn get_mut(&mut self, id: Id<T>) -> Option<&mut T> {
        let i = id.index() as usize;
        if i >= self.items.len() || self.removed[i] {
            return None;
        }
        self.revisions[i].own = self.revisions[i].own.wrapping_add(1);
        self.contents[i] = self.contents[i].wrapping_add(1);
        Some(&mut self.items[i])
    }

    /// [`Table::get_mut`] for an edit of what is kept with the entry but not of what it is
    /// (the link to its mesh or polyline): an edit all the same (its own revision moves on),
    /// but values computed from what it is stay good.
    pub(crate) fn get_mut_keeping_content(&mut self, id: Id<T>) -> Option<&mut T> {
        let i = id.index() as usize;
        if i >= self.items.len() || self.removed[i] {
            return None;
        }
        self.revisions[i].own = self.revisions[i].own.wrapping_add(1);
        Some(&mut self.items[i])
    }

    /// How often what the entry is was edited (also available for removed entries).
    pub(crate) fn content(&self, id: Id<T>) -> Option<u32> {
        self.contents.get(id.index() as usize).copied()
    }

    /// Persistent counter of the entry (also available for removed entries).
    pub fn stable_id(&self, id: Id<T>) -> Option<u32> {
        self.stable_ids.get(id.index() as usize).copied()
    }

    /// Change counters of the entry (also available for removed entries).
    pub fn revision(&self, id: Id<T>) -> Option<Revision> {
        self.revisions.get(id.index() as usize).copied()
    }

    pub fn revision_mut(&mut self, id: Id<T>) -> Option<&mut Revision> {
        self.revisions.get_mut(id.index() as usize)
    }

    /// Soft-removes the entry; returns false if it was already removed or out of range.
    pub fn remove(&mut self, id: Id<T>) -> bool {
        let i = id.index() as usize;
        if i >= self.items.len() || self.removed[i] {
            return false;
        }
        self.removed[i] = true;
        self.active -= 1;
        self.revisions[i].own = self.revisions[i].own.wrapping_add(1);
        self.contents[i] = self.contents[i].wrapping_add(1);
        true
    }

    /// [`Table::remove`], then `release` on the removed entry to free what it holds, and the
    /// slot is reused by a later [`Table::push`]: for entries nothing live refers to.
    pub(crate) fn remove_releasing(&mut self, id: Id<T>, release: impl FnOnce(&mut T)) -> bool {
        if !self.remove(id) {
            return false;
        }
        let i = id.index() as usize;
        release(&mut self.items[i]);
        self.released[i] = true;
        self.free.push(id.index());
        true
    }

    /// Overwrites the entry, reviving it if it was removed. Panics if `id` is out of range.
    pub fn replace(&mut self, id: Id<T>, item: T) {
        let i = id.index() as usize;
        assert!(i < self.items.len(), "replace of out-of-range {id:?}");
        if self.removed[i] {
            self.removed[i] = false;
            self.released[i] = false;
            self.active += 1;
        }
        self.items[i] = item;
        self.revisions[i].own = self.revisions[i].own.wrapping_add(1);
        self.contents[i] = self.contents[i].wrapping_add(1);
    }

    pub fn is_removed(&self, id: Id<T>) -> bool {
        self.removed.get(id.index() as usize).copied().unwrap_or(true)
    }

    /// Number of slots, removed ones included (the exclusive upper bound of indices).
    pub fn capacity_len(&self) -> usize {
        self.items.len()
    }

    /// Number of entries that are not removed.
    pub fn len(&self) -> usize {
        self.active
    }

    pub fn is_empty(&self) -> bool {
        self.active == 0
    }

    /// Active entries in index order.
    pub fn iter(&self) -> impl Iterator<Item = (Id<T>, &T)> + '_ {
        self.items
            .iter()
            .enumerate()
            .filter(|(i, _)| !self.removed[*i])
            .map(|(i, t)| (Id::from_index(i as u32), t))
    }

    /// All entries in index order, with whether each is removed.
    pub fn iter_all(&self) -> impl Iterator<Item = (Id<T>, &T, bool)> + '_ {
        self.items
            .iter()
            .enumerate()
            .map(|(i, t)| (Id::from_index(i as u32), t, self.removed[i]))
    }

    /// Every entry for rewriting ids stored in it (not an edit: revisions stay).
    pub(crate) fn items_mut(&mut self) -> impl Iterator<Item = &mut T> + '_ {
        self.items.iter_mut()
    }

    /// The live entry whose stable id is `stable_id`.
    pub(crate) fn find_stable(&self, stable_id: u32) -> Option<Id<T>> {
        let i = self.stable_ids.iter().position(|&s| s == stable_id)?;
        (!self.removed[i]).then(|| Id::from_index(i as u32))
    }

    /// Drops removed slots and renumbers the rest, preserving order. Returns the old-index ->
    /// new-id map (`Id::INVALID` for removed entries); ids stored elsewhere must be remapped
    /// with it.
    pub fn compact(&mut self) -> Vec<Id<T>> {
        let mut remap = Vec::with_capacity(self.items.len());
        let mut next = 0u32;
        for &removed in &self.removed {
            if removed {
                remap.push(Id::INVALID);
            } else {
                remap.push(Id::from_index(next));
                next += 1;
            }
        }
        let mut i = 0;
        let removed = std::mem::take(&mut self.removed);
        self.items.retain(|_| {
            let keep = !removed[i];
            i += 1;
            keep
        });
        let mut i = 0;
        self.stable_ids.retain(|_| {
            let keep = !removed[i];
            i += 1;
            keep
        });
        let mut i = 0;
        self.revisions.retain(|_| {
            let keep = !removed[i];
            i += 1;
            keep
        });
        let mut i = 0;
        self.contents.retain(|_| {
            let keep = !removed[i];
            i += 1;
            keep
        });
        self.removed = vec![false; self.items.len()];
        self.released = vec![false; self.items.len()];
        self.free.clear();
        remap
    }
}

impl<T> Index<Id<T>> for Table<T> {
    type Output = T;
    /// Panics on an invalid or removed id.
    fn index(&self, id: Id<T>) -> &T {
        self.get(id).unwrap_or_else(|| panic!("invalid or removed {id:?}"))
    }
}

impl<T> IndexMut<Id<T>> for Table<T> {
    fn index_mut(&mut self, id: Id<T>) -> &mut T {
        self.get_mut(id).unwrap_or_else(|| panic!("invalid or removed {id:?}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn released_slots_are_reused_unless_revived() {
        let mut t = Table::default();
        let a = t.push("a", 0);
        let b = t.push("b", 1);
        let c = t.push("c", 2);
        assert!(t.remove_releasing(a, |v| *v = ""));
        assert!(t.remove_releasing(b, |v| *v = ""));
        t.replace(a, "a again");
        // `b` is reused (with a new stable id and counters moved on); revived `a` is not.
        let content = t.content(b);
        assert_eq!(t.push("d", 3), b);
        assert_eq!((t[b], t.stable_id(b)), ("d", Some(3)));
        assert_ne!(t.content(b), content);
        assert_eq!(t.push("e", 4), Id::from_index(3));
        assert_eq!(t.iter().map(|(_, v)| *v).collect::<Vec<_>>(), ["a again", "d", "c", "e"]);
        // A plain removal is not reused.
        assert!(t.remove(c));
        assert_eq!(t.push("f", 5), Id::from_index(4));
    }

    #[test]
    fn soft_removal_keeps_other_ids_valid() {
        let mut t = Table::default();
        let a = t.push("a", 0);
        let b = t.push("b", 1);
        let c = t.push("c", 2);
        assert!(t.remove(b));
        assert!(!t.remove(b));
        assert_eq!(t.len(), 2);
        assert_eq!(t[a], "a");
        assert_eq!(t[c], "c");
        assert!(t.get(b).is_none());
        assert_eq!(t.stable_id(b), Some(1));
        assert_eq!(t.iter().map(|(_, v)| *v).collect::<Vec<_>>(), ["a", "c"]);
    }

    #[test]
    fn compact_renumbers_in_order() {
        let mut t = Table::default();
        let ids: Vec<_> = (0..4).map(|i| t.push(i * 10, i as u32)).collect();
        t.remove(ids[1]);
        let remap = t.compact();
        assert_eq!(remap[1], Id::INVALID);
        assert_eq!(remap[3], Id::from_index(2));
        assert_eq!(t[remap[3]], 30);
        assert_eq!(t.stable_id(remap[3]), Some(3));
        assert_eq!(t.capacity_len(), 3);
    }

    #[test]
    fn replace_revives_removed_entry_in_place() {
        let mut t = Table::default();
        let a = t.push(1, 0);
        t.remove(a);
        assert_eq!(t.len(), 0);
        t.replace(a, 2);
        assert_eq!(t.len(), 1);
        assert_eq!(t[a], 2);
        assert_eq!(t.stable_id(a), Some(0));
    }

    #[test]
    #[should_panic(expected = "invalid or removed")]
    fn indexing_removed_entry_panics() {
        let mut t = Table::default();
        let a = t.push(1, 0);
        t.remove(a);
        let _ = t[a];
    }
}
