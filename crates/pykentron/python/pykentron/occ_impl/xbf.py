"""OCCT's binary OCAF documents (XBF, "BinXCAF", format version 12, as OCCT 7.6 and later
write them), as CadQuery writes and reads assemblies in them.

The file: "BINFILE", a header of 13 integers, the information section (format, date, the
attribute types used), comments, a table of sections, then the label tree: each label its
tag, its size, its attributes (type, id, length, payload: values aligned to 4 bytes) and its
child labels. Shapes and placements are written in the attributes holding them, as OCCT's
shape streams (kentron's `ShapeStreams`: each shape where first met, later uses referring
back to it). Everything little-endian. Written from the same label tree as the XML format's
(`ocaf._build`), read into the same questions CadQuery's `_importDoc` asks."""

from __future__ import annotations

import datetime
import struct
import uuid

from .. import _kentron as native
from ..shapes import Shape
from . import ocaf

MAGIC = b"BINFILE"
END_ATTRIBUTES = -1
END_LABEL = -2
NAME_GUID = ocaf._NAME

# Attribute types whose payload is empty (the tools).
EMPTY_TYPES = {
    "XCAFDoc_DocumentTool", "XCAFDoc_ShapeTool", "XCAFDoc_ColorTool", "XCAFDoc_LayerTool",
    "XCAFDoc_DimTolTool", "XCAFDoc_MaterialTool", "XCAFDoc_ViewTool", "XCAFDoc_ClippingPlaneTool",
    "XCAFDoc_NotesTool", "XCAFDoc_VisMaterialTool",
}
# Attributes whose shapes follow their payload as a stream (its size first).
STREAM_TYPES = {"TNaming_NamedShape", "XCAFDoc_Location"}


def _guid_bytes(s: str) -> bytes:
    return uuid.UUID(s).bytes_le


def _guid(b: bytes) -> str:
    return str(uuid.UUID(bytes_le=bytes(b)))


class _Payload:
    """An attribute's payload (OCCT's BinObjMgt_Persistent): values aligned to 4 bytes but
    characters."""

    def __init__(self, data: bytes = b""):
        self.b = bytearray(data)
        self.at = 0

    # Writing.
    def put_int(self, v: int):
        while len(self.b) % 4:
            self.b.append(0)
        self.b += struct.pack("<i", v)

    def put_real(self, v: float):
        while len(self.b) % 4:
            self.b.append(0)
        self.b += struct.pack("<d", v)

    def put_shortreal(self, v: float):
        while len(self.b) % 4:
            self.b.append(0)
        self.b += struct.pack("<f", v)

    def put_char(self, c: str):
        self.b += c.encode("latin-1")

    def put_guid(self, g: str):
        while len(self.b) % 4:
            self.b.append(0)
        self.b += _guid_bytes(g)

    def put_extstring(self, s: str):
        while len(self.b) % 4:
            self.b.append(0)
        self.b += s.encode("utf-16-le") + b"\0\0"

    # Reading.
    def _skip_to_4(self):
        if self.at % 4:
            self.at += 4 - self.at % 4

    def int(self) -> int:
        self._skip_to_4()
        v = struct.unpack_from("<i", self.b, self.at)[0]
        self.at += 4
        return v

    def real(self) -> float:
        self._skip_to_4()
        v = struct.unpack_from("<d", self.b, self.at)[0]
        self.at += 8
        return v

    def shortreal(self) -> float:
        self._skip_to_4()
        v = struct.unpack_from("<f", self.b, self.at)[0]
        self.at += 4
        return v

    def char(self) -> str:
        v = self.b[self.at:self.at + 1].decode("latin-1")
        self.at += 1
        return v

    def guid(self) -> str:
        self._skip_to_4()
        v = _guid(self.b[self.at:self.at + 16])
        self.at += 16
        return v

    def extstring(self) -> str:
        self._skip_to_4()
        end = self.at
        while self.b[end:end + 2] != b"\0\0":
            end += 2
        s = self.b[self.at:end].decode("utf-16-le")
        self.at = end + 2
        return s


def _str(s: str) -> bytes:
    b = s.encode("latin-1")
    return struct.pack("<i", len(b)) + b


def _extstr(s: str) -> bytes:
    b = s.encode("utf-16-le")
    return struct.pack("<i", len(b) // 2) + b


def _to_linear(c: float) -> float:
    """An sRGB component as the linear one OCCT keeps (single precision)."""
    return struct.unpack("f", struct.pack("f", c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4))[0]


# Writing.


def _entry(label) -> str:
    return label.entry()


def _attributes(doc: "ocaf._Document"):
    """The document's labels as `(tag, [attribute dicts], children)`, attributes as OCCT
    writes them."""

    def tree_links(node: "ocaf._Node"):
        father = node.fathers[0] if node.fathers else None
        siblings = father.children if father is not None else []
        k = siblings.index(node) if father is not None else -1
        return {
            "father": _entry(father.label) if father is not None else None,
            "next": _entry(siblings[k + 1].label) if father is not None and k + 1 < len(siblings) else None,
            "previous": _entry(siblings[k - 1].label) if father is not None and k > 0 else None,
            "first": _entry(node.children[0].label) if node.children else None,
        }

    def convert(label):
        attributes = []
        for kind, value in label.attributes:
            key = kind[0] if isinstance(kind, tuple) else kind
            if key == "raw":
                attributes.append({"type": value})
            elif key == "name":
                attributes.append({"type": "TDataStd_Name", "value": value})
            elif key == "uattr":
                attributes.append({"type": "TDataStd_UAttribute", "guid": value})
            elif key == "color":
                rgba = label.rgba
                attributes.append({"type": "XCAFDoc_Color", "rgb": [_to_linear(c) for c in rgba[:3]], "alpha": rgba[3]})
            elif key == "location":
                attributes.append({"type": "XCAFDoc_Location", "location": value})
            elif key == "shape":
                attributes.append({"type": "TNaming_NamedShape", "version": 1 if label.parent is doc.shapes else 0, "shape": value})
            elif key == "tree":
                attributes.append({"type": "TDataStd_TreeNode", "treeid": value.tree, "links": tree_links(value), "_node": value})
            elif key == "graph":
                attributes.append({"type": "XCAFDoc_GraphNode", "graphid": value.tree, "fathers": [_entry(f.label) for f in value.fathers], "children": [_entry(c.label) for c in value.children], "_node": value})
        if label.next_tag:
            attributes.append({"type": "TDF_TagSource", "value": label.next_tag})
        return {"tag": label.tag, "attributes": attributes, "children": [convert(c) for c in label.children.values()], "_label": label}

    return convert(doc.root)


def _walk(label, path=()):
    path = path + (label["tag"],)
    yield path, label
    for c in label["children"]:
        yield from _walk(c, path)


def _ids(root):
    """OCCT's attribute ids (its relocation table): each attribute numbered as written,
    the tree and graph nodes it links to numbered where first named."""
    by_entry = {":".join(map(str, path)): label["attributes"] for path, label in _walk(root)}
    table = {}

    def add(a):
        if id(a) not in table:
            table[id(a)] = len(table) + 1
        return table[id(a)]

    def find(entry, type_name, key, value):
        return next(a for a in by_entry[entry] if a["type"] == type_name and a.get(key) == value)

    result = {}
    for _, label in _walk(root):
        for a in label["attributes"]:
            own = add(a)
            extra = None
            if a["type"] == "TDataStd_TreeNode":
                extra = [-1 if a["links"][k] is None else add(find(a["links"][k], "TDataStd_TreeNode", "treeid", a["treeid"])) for k in ("father", "next", "previous", "first")]
            elif a["type"] == "XCAFDoc_GraphNode":
                extra = ([add(find(e, "XCAFDoc_GraphNode", "graphid", a["graphid"])) for e in a["fathers"]], [add(find(e, "XCAFDoc_GraphNode", "graphid", a["graphid"])) for e in a["children"]])
            result[id(a)] = (own, extra)
    return result


def _payload(a, extra) -> bytes:
    w = _Payload()
    t = a["type"]
    if t == "TDataStd_Name":
        w.put_extstring(a["value"])
        w.put_guid(NAME_GUID)
    elif t == "TDataStd_TreeNode":
        for v in extra:
            w.put_int(v)
        w.put_guid(a["treeid"])
    elif t == "TDataStd_UAttribute":
        w.put_guid(a["guid"])
    elif t == "TDF_TagSource":
        w.put_int(a["value"])
    elif t == "XCAFDoc_Color":
        for v in a["rgb"]:
            w.put_real(v)
        w.put_shortreal(a["alpha"])
    elif t == "XCAFDoc_GraphNode":
        fathers, children = extra
        for v in fathers:
            w.put_int(v)
        w.put_int(-1)
        for v in children:
            w.put_int(v)
        w.put_int(-1)
        w.put_guid(a["graphid"])
    elif t == "TNaming_NamedShape":
        # (One shape, made from nothing: "primitive".)
        w.put_int(1)
        w.put_int(a["version"])
        w.put_char("P")
    elif t in EMPTY_TYPES or t == "XCAFDoc_Location":
        pass
    else:
        raise ValueError(f"cannot write attribute type {t}")
    return bytes(w.b)


def write_xbf(assy, path: str) -> bool:
    doc = ocaf._build(assy)
    # (Colour labels keep their components for the binary format, which stores them.)
    for label, rgba in doc.colors:
        label.rgba = rgba
    root = _attributes(doc)
    types = []
    for _, label in _walk(root):
        for a in label["attributes"]:
            if a["type"] not in types:
                types.append(a["type"])
    user_info = ["FILE_FORMAT: BinXCAF", "REFERENCE_COUNTER: 0", "MODIFICATION_COUNTER: 1", "START_TYPES"] + types + ["END_TYPES"]
    info = (struct.pack("<i", 1) + _str("12") + _str(datetime.date.today().strftime("%m/%d/%Y")) + _str("") + _str("1")
            + _extstr("") + _str("") + _extstr("") + struct.pack("<i", len(user_info)) + b"".join(_str(s) for s in user_info))
    comments = struct.pack("<i", 0)
    binfo = len(MAGIC)
    einfo = binfo + 13 * 4 + len(info)
    ecomment = einfo + len(comments)
    header = [0x04030201, binfo, einfo, einfo, ecomment, -1, -1, -1, -1, -1, -1, -1, ecomment]
    out = bytearray(MAGIC + struct.pack("<13i", *header) + info + comments)
    # The table of sections: only its end (":").
    out += struct.pack("<i", 4) + b":\0\0\0" + struct.pack("<3Q", 0, 0, 0)
    ids = _ids(root)
    type_ids = {t: i + 1 for i, t in enumerate(types)}
    streams = native.ShapeStreams()

    def write(label):
        out.extend(struct.pack("<i", label["tag"]))
        size_at = len(out)
        out.extend(b"\0" * 8)
        for a in label["attributes"]:
            own, extra = ids[id(a)]
            payload = _payload(a, extra)
            if a["type"] in STREAM_TYPES:
                # (A negative id: a stream follows, its size first.)
                out.extend(struct.pack("<iii", type_ids[a["type"]], -own, len(payload)) + payload)
                at = len(out)
                out.extend(b"\0" * 8)
                if a["type"] == "XCAFDoc_Location":
                    loc = a["location"]
                    rows = loc.rows if loc is not None else [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]
                    out.extend(streams.location(rows, len(out)))
                else:
                    out.extend(streams.shape(a["shape"].wrapped, len(out)))
                struct.pack_into("<Q", out, at, len(out) - at)
            else:
                out.extend(struct.pack("<iii", type_ids[a["type"]], own, len(payload)) + payload)
        out.extend(struct.pack("<i", END_ATTRIBUTES))
        for c in label["children"]:
            write(c)
        out.extend(struct.pack("<i", END_LABEL))
        struct.pack_into("<Q", out, size_at, len(out) - size_at)

    write(root)
    with open(path, "wb") as f:
        f.write(bytes(out))
    return True


# Reading.


class _ReadBin(ocaf._Read):
    """A binary document read: its labels, attributes and shapes, asked as an XML one's."""

    def __init__(self, path: str):
        with open(path, "rb") as f:
            data = f.read()
        if data[:7] != MAGIC:
            raise ValueError(f"not an OCAF binary document: {path}")
        at = 7 + 13 * 4
        n, = struct.unpack_from("<i", data, at)
        at += 4

        def read_str():
            nonlocal at
            k, = struct.unpack_from("<i", data, at)
            at += 4
            s = data[at:at + max(k, 0)].decode("latin-1")
            at += max(k, 0)
            return s

        def read_extstr():
            nonlocal at
            k, = struct.unpack_from("<i", data, at)
            at += 4
            s = data[at:at + 2 * max(k, 0)].decode("utf-16-le")
            at += 2 * max(k, 0)
            return s

        version = read_str()
        read_str(), read_str(), read_str()
        read_extstr(), read_str(), read_extstr()
        count, = struct.unpack_from("<i", data, at)
        at += 4
        user_info = [read_str() for _ in range(count)]
        comments, = struct.unpack_from("<i", data, at)
        at += 4
        for _ in range(comments):
            read_extstr()
        if not version.isdigit() or int(version) < 12:
            raise NotImplementedError(f"binary OCAF documents of format version {version} (12 and later are read)")
        types = user_info[user_info.index("START_TYPES") + 1:user_info.index("END_TYPES")]
        while True:
            k, = struct.unpack_from("<i", data, at)
            at += 4
            name = data[at:at + k].rstrip(b"\0").decode("latin-1")
            at += k + 24
            if name == ":":
                break
        self.labels = {}
        self.nodes = {}
        self.text = ""
        self._table = []
        streams_at = []  # (entry, kind, position)
        by_id = {}  # attribute id -> (entry, kind, attribute)

        def label(parent):
            nonlocal at
            tag, = struct.unpack_from("<i", data, at)
            at += 4 + 8
            entry = f"{parent}:{tag}" if parent is not None else str(tag)
            record = {"entry": entry, "tag": tag, "parent": parent, "children": [], "name": None, "shape": None, "location": None, "uattr": set(), "trees": {}, "graphs": {}, "color": None, "rgba": None}
            self.labels[entry] = record
            if parent is not None:
                self.labels[parent]["children"].append(entry)
            while True:
                t, = struct.unpack_from("<i", data, at)
                at += 4
                if t == END_ATTRIBUTES:
                    break
                own, length = struct.unpack_from("<ii", data, at)
                at += 8
                p = _Payload(data[at:at + length])
                at += length
                kind = types[t - 1]
                if own < 0:
                    size, = struct.unpack_from("<Q", data, at)
                    streams_at.append((entry, kind, at + 8))
                    if kind == "TNaming_NamedShape" and (p.int(), p.int(), p.char())[0] != 1:
                        raise NotImplementedError("named shapes of other than one shape")
                    at += size
                elif kind == "TDataStd_Name":
                    record["name"] = p.extstring()
                elif kind == "TDataStd_UAttribute":
                    record["uattr"].add(p.guid())
                elif kind == "XCAFDoc_Color":
                    record["rgba"] = (p.real(), p.real(), p.real(), p.shortreal())
                elif kind == "TDataStd_TreeNode":
                    links = [p.int() for _ in range(4)]
                    node = {"id": abs(own), "label": entry, "tree": p.guid(), "links": links, "father": None, "children": [], "fathers": []}
                    record["trees"][node["tree"]] = node
                    by_id[abs(own)] = node
                elif kind == "XCAFDoc_GraphNode":
                    fathers = []
                    while (v := p.int()) != -1:
                        fathers.append(v)
                    children = []
                    while (v := p.int()) != -1:
                        children.append(v)
                    node = {"id": abs(own), "label": entry, "tree": p.guid(), "fathers_ids": fathers, "children": children, "fathers": []}
                    record["graphs"][node["tree"]] = node
                    by_id[abs(own)] = node
            while struct.unpack_from("<i", data, at)[0] != END_LABEL:
                label(entry)
            at += 4

        label(None)
        for node in by_id.values():
            if "links" in node and node["links"][0] >= 0:
                node["father"] = by_id.get(node["links"][0])
            if "fathers_ids" in node:
                node["fathers"] = [str(i) for i in node["fathers_ids"] if i in by_id]
        self.nodes = {str(i): n for i, n in by_id.items()}
        shape_at = [p for _, kind, p in streams_at if kind == "TNaming_NamedShape"]
        location_at = [p for _, kind, p in streams_at if kind == "XCAFDoc_Location"]
        shapes, rows = native.read_shape_streams(data, shape_at, location_at)
        self.shapes_cache = {}
        for (entry, _, _), handle in zip([s for s in streams_at if s[1] == "TNaming_NamedShape"], shapes):
            self.labels[entry]["shape"] = True
            if handle is not None:
                self.shapes_cache[entry] = Shape.cast(handle)
        for (entry, _, _), r in zip([s for s in streams_at if s[1] == "XCAFDoc_Location"], rows):
            self._table.append([list(x) for x in r])
            self.labels[entry]["location"] = len(self._table)

    def color(self, entry: str, trees=(ocaf._COLOR_SURF, ocaf._COLOR_GEN)):
        from ..geom import Color

        for tree in trees:
            father = self.father(entry, tree)
            if father is not None and self.labels[father]["rgba"] is not None:
                r, g, b, a = self.labels[father]["rgba"]
                return Color(r, g, b, a, False)
        return None

    def shapes(self, entries: list) -> list:
        return [self.shapes_cache.get(e) for e in entries]


def read_xbf(assy, path: str) -> None:
    """Fills `assy` from a binary OCAF document, as CadQuery's `_importDoc` does."""
    ocaf.read_document(assy, _ReadBin(path), path)
