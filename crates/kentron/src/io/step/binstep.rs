//! `.binstep`: a STEP file's exchange structure (ISO 10303-21) in binary, value for value
//! the same (header records, instances by number and in order, records of complex
//! instances, every parameter), for any schema; reading it builds the same [`Exchange`] as
//! parsing the text, without the text's tokenizing and number conversion.
//!
//! Layout (integers as LEB128 varints, signed ones zigzagged):
//!
//! ```text
//! "BINSTEP\x01"
//! names:     count, then each: length, UTF-8 bytes    (type and enumeration names)
//! header:    count, then each record
//! instances: count, then each: number (less the previous one's, signed),
//!            record count, records
//! record:    name (index into names), list
//! list:      count, then each value
//! value:     tag byte, then by tag:
//!            0 $   1 *   2 integer (signed)   3 real (8 bytes, little-endian IEEE 754)
//!            4 string (length, UTF-8)   5 binary (length, hexadecimal digits)
//!            6 enumeration (name)   7 reference (number less the instance's own, signed)
//!            8 list   9 typed value (name, list)
//! ```

use crate::error::{Error, Result};
use crate::hashing::FastMap as HashMap;

use super::parse::{Exchange, Name, Stretch, Value};

/// The first bytes of a `.binstep` file (the last one its version).
pub const MAGIC: &[u8; 8] = b"BINSTEP\x01";

/// True if `bytes` start as a `.binstep` file does.
pub fn is_binstep(bytes: &[u8]) -> bool {
    bytes.starts_with(MAGIC)
}

/// The exchange as `.binstep` bytes.
pub fn encode(x: &Exchange) -> Vec<u8> {
    let mut e = Encoder { out: Vec::with_capacity(x.values.len() * 4), names: HashMap::default(), table: Vec::new() };
    // The names first (in order of first use): records', enumerations' and typed values'.
    for (name, params) in x.header.iter().chain(&x.records) {
        e.name(name);
        e.collect(x, *params);
    }
    write_varint(&mut e.out, x.header.len() as u64);
    for (name, params) in &x.header {
        e.record(x, name, *params, 0);
    }
    write_varint(&mut e.out, x.instances.len() as u64);
    let mut previous = 0i64;
    for &(id, records) in &x.instances {
        write_varint(&mut e.out, zigzag(id as i64 - previous));
        previous = id as i64;
        let records = &x.records[records.range()];
        write_varint(&mut e.out, records.len() as u64);
        for (name, params) in records {
            e.record(x, name, *params, id as i64);
        }
    }
    let mut out = Vec::with_capacity(e.out.len() + 32 * e.table.len() + 16);
    out.extend_from_slice(MAGIC);
    write_varint(&mut out, e.table.len() as u64);
    for name in &e.table {
        write_varint(&mut out, name.len() as u64);
        out.extend_from_slice(name.as_bytes());
    }
    out.extend_from_slice(&e.out);
    out
}

struct Encoder<'a> {
    out: Vec<u8>,
    names: HashMap<&'a str, u32>,
    table: Vec<&'a str>,
}

impl<'a> Encoder<'a> {
    /// The name's index, given one at its first use.
    fn name(&mut self, name: &'a str) -> u32 {
        let next = self.table.len() as u32;
        *self.names.entry(name).or_insert_with(|| {
            self.table.push(name);
            next
        })
    }

    /// The names used inside a list (enumerations, typed values), in order.
    fn collect(&mut self, x: &'a Exchange, list: Stretch) {
        for v in &x.values[list.range()] {
            match v {
                Value::Enumeration(e) => {
                    self.name(e);
                }
                Value::List(l) => self.collect(x, *l),
                Value::Typed(t) => {
                    self.name(&t.0);
                    self.collect(x, t.1);
                }
                _ => {}
            }
        }
    }

    fn record(&mut self, x: &'a Exchange, name: &'a str, params: Stretch, id: i64) {
        let k = self.name(name);
        write_varint(&mut self.out, k as u64);
        self.list(x, params, id);
    }

    fn list(&mut self, x: &'a Exchange, list: Stretch, id: i64) {
        write_varint(&mut self.out, list.len as u64);
        for v in &x.values[list.range()] {
            self.value(x, v, id);
        }
    }

    fn value(&mut self, x: &'a Exchange, v: &'a Value, id: i64) {
        match v {
            Value::Unset => self.out.push(0),
            Value::Derived => self.out.push(1),
            Value::Integer(i) => {
                self.out.push(2);
                write_varint(&mut self.out, zigzag(*i));
            }
            Value::Real(r) => {
                self.out.push(3);
                self.out.extend_from_slice(&r.to_le_bytes());
            }
            Value::String(s) => {
                self.out.push(4);
                write_varint(&mut self.out, s.len() as u64);
                self.out.extend_from_slice(s.as_bytes());
            }
            Value::Binary(b) => {
                self.out.push(5);
                write_varint(&mut self.out, b.len() as u64);
                self.out.extend_from_slice(b.as_bytes());
            }
            Value::Enumeration(e) => {
                self.out.push(6);
                let k = self.name(e);
                write_varint(&mut self.out, k as u64);
            }
            Value::Reference(r) => {
                self.out.push(7);
                write_varint(&mut self.out, zigzag(*r as i64 - id));
            }
            Value::List(l) => {
                self.out.push(8);
                self.list(x, *l, id);
            }
            Value::Typed(t) => {
                self.out.push(9);
                self.record(x, &t.0, t.1, id);
            }
        }
    }
}

/// The exchange in `.binstep` bytes.
pub fn decode(bytes: &[u8]) -> Result<Exchange> {
    if !is_binstep(bytes) {
        return Err(invalid("not a .binstep file"));
    }
    if bytes.len() >= u32::MAX as usize {
        return Err(invalid("files of 4 GB and more are not read"));
    }
    let mut d = Decoder { s: bytes, i: MAGIC.len(), names: Vec::new(), values: Vec::new() };
    let count = d.count()?;
    for _ in 0..count {
        let name = Name::from(d.text()?);
        d.names.push(name);
    }
    let mut x = Exchange::default();
    // (About a value per 5 bytes.)
    d.values.reserve(bytes.len() / 5);
    let count = d.count()?;
    for _ in 0..count {
        let record = d.record(0)?;
        x.header.push(record);
    }
    let count = d.count()?;
    x.instances.reserve(count);
    x.records.reserve(count);
    x.places.reserve(count);
    let mut by_name: Vec<Vec<u32>> = vec![Vec::new(); d.names.len()];
    let mut id = 0i64;
    for _ in 0..count {
        id += unzigzag(d.varint()?);
        if id < 0 {
            return Err(invalid("a negative instance number"));
        }
        let place = x.instances.len() as u32;
        let first = x.records.len() as u32;
        let records = d.count()?;
        for _ in 0..records {
            let k = d.count()?;
            let name = d.names.get(k).cloned().ok_or_else(|| invalid("a name out of range"))?;
            let list = &mut by_name[k];
            if list.last() != Some(&place) {
                list.push(place);
            }
            let params = d.list(id)?;
            x.records.push((name, params));
        }
        // (A number given twice: the last instance, as when parsing.)
        x.places.insert(id as usize, place);
        x.instances.push((id as usize, Stretch { start: first, len: x.records.len() as u32 - first }));
    }
    if d.i != bytes.len() {
        return Err(invalid("bytes after the instances"));
    }
    x.values = std::mem::take(&mut d.values);
    for (k, places) in by_name.into_iter().enumerate() {
        if !places.is_empty() {
            x.types.insert(d.names[k].clone(), places);
        }
    }
    Ok(x)
}

struct Decoder<'a> {
    s: &'a [u8],
    i: usize,
    names: Vec<Name>,
    /// The lists' values: each list's slots taken when it starts, then filled (lists within
    /// it taking theirs after).
    values: Vec<Value>,
}

impl Decoder<'_> {
    fn take(&mut self, n: usize) -> Result<&[u8]> {
        let end = self.i.checked_add(n).filter(|&e| e <= self.s.len()).ok_or_else(|| invalid("a file cut short"))?;
        let out = &self.s[self.i..end];
        self.i = end;
        Ok(out)
    }

    fn byte(&mut self) -> Result<u8> {
        let b = *self.s.get(self.i).ok_or_else(|| invalid("a file cut short"))?;
        self.i += 1;
        Ok(b)
    }

    #[inline]
    fn varint(&mut self) -> Result<u64> {
        // (Mostly one byte.)
        if let Some(&b) = self.s.get(self.i) {
            if b < 0x80 {
                self.i += 1;
                return Ok(u64::from(b));
            }
        }
        let mut out = 0u64;
        for shift in (0..64).step_by(7) {
            let b = self.byte()?;
            out |= u64::from(b & 0x7f) << shift;
            if b & 0x80 == 0 {
                return Ok(out);
            }
        }
        Err(invalid("a number too long"))
    }

    /// A count or index (bounded by the file's size, so that a bad one fails rather than
    /// allocating).
    fn count(&mut self) -> Result<usize> {
        let n = self.varint()?;
        if n > self.s.len() as u64 {
            return Err(invalid("a count beyond the file"));
        }
        Ok(n as usize)
    }

    fn name(&mut self) -> Result<Name> {
        let k = self.count()?;
        self.names.get(k).cloned().ok_or_else(|| invalid("a name out of range"))
    }

    fn text(&mut self) -> Result<Box<str>> {
        let n = self.count()?;
        let bytes = self.take(n)?;
        Ok(std::str::from_utf8(bytes).map_err(|_| invalid("text that is not UTF-8"))?.into())
    }

    fn record(&mut self, id: i64) -> Result<(Name, Stretch)> {
        let name = self.name()?;
        Ok((name, self.list(id)?))
    }

    fn list(&mut self, id: i64) -> Result<Stretch> {
        let n = self.count()?;
        let start = self.values.len();
        self.values.resize(start + n, Value::Unset);
        for k in 0..n {
            let v = self.value(id)?;
            self.values[start + k] = v;
        }
        Ok(Stretch { start: start as u32, len: n as u32 })
    }

    fn value(&mut self, id: i64) -> Result<Value> {
        Ok(match self.byte()? {
            0 => Value::Unset,
            1 => Value::Derived,
            2 => Value::Integer(unzigzag(self.varint()?)),
            3 => Value::Real(f64::from_le_bytes(self.take(8)?.try_into().expect("8 bytes"))),
            4 => Value::String(self.text()?),
            5 => Value::Binary(self.text()?),
            6 => Value::Enumeration(self.name()?),
            7 => {
                let r = id + unzigzag(self.varint()?);
                if r < 0 {
                    return Err(invalid("a negative reference"));
                }
                Value::Reference(r as usize)
            }
            8 => Value::List(self.list(id)?),
            9 => Value::Typed(Box::new(self.record(id)?)),
            tag => return Err(invalid(&format!("an unknown value tag {tag}"))),
        })
    }
}

fn invalid(what: &str) -> Error {
    Error::InvalidInput(format!(".binstep: {what}"))
}

fn write_varint(out: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        out.push((v as u8) | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

fn zigzag(v: i64) -> u64 {
    ((v << 1) ^ (v >> 63)) as u64
}

fn unzigzag(v: u64) -> i64 {
    ((v >> 1) as i64) ^ -((v & 1) as i64)
}

#[cfg(test)]
mod tests {
    use super::super::parse::parse;
    use super::*;

    const SAMPLE: &[u8] = b"ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('a test'),'2;1');\nFILE_NAME('caf\\X2\\00E9\\X0\\','2026-10-09',('x'),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n#1 = CARTESIAN_POINT('it''s',(0.,1.5,-2.E-03));\n#2 = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );\n#3 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#2,'distance_accuracy_value','');\n#10 = DIRECTION('',(0.,0.,1.)) /* a comment */;\n#7 = SOMETHING($,\"0F\",-12,(#1,(#10,#3)),.T.);\nENDSEC;\nEND-ISO-10303-21;\n";

    #[test]
    fn round_trips_through_binary_and_text() {
        let x = parse(SAMPLE).unwrap();
        let text = x.to_text();
        let bytes = encode(&x);
        assert!(is_binstep(&bytes));
        let y = decode(&bytes).unwrap();
        assert_eq!(y.to_text(), text);
        // The text written reads back the same, and so does its binary.
        let z = parse(text.as_bytes()).unwrap();
        assert_eq!(z.to_text(), text);
        assert_eq!(encode(&z), bytes);
        assert_eq!(y.header.len(), 3);
        assert_eq!(y.ids().collect::<Vec<_>>(), vec![1, 2, 3, 10, 7]);
        assert_eq!(y.get(7).unwrap().params()[1], Value::Binary("0F".into()));
        assert_eq!(y.of_type("SI_UNIT").map(|(id, _)| id).collect::<Vec<_>>(), vec![2]);
    }

    #[test]
    fn cut_short_or_garbled_files_fail() {
        let bytes = encode(&parse(SAMPLE).unwrap());
        for n in [8, 20, bytes.len() / 2, bytes.len() - 1] {
            assert!(decode(&bytes[..n]).is_err(), "cut at {n}");
        }
        assert!(decode(b"BINSTEP\x01\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\x01").is_err());
    }
}
