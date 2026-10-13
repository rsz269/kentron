//! The exchange-structure syntax of STEP files (ISO 10303-21): the header's records and the
//! data section's instances parsed into values, by instance number; and written back.

use crate::hashing::FastMap as HashMap;

use crate::error::{Error, Result};

/// A parameter value (kept small: files hold millions).
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Integer(i64),
    Real(f64),
    String(Box<str>),
    /// `"..."`: binary, as its hexadecimal digits (the first telling the unused bits).
    Binary(Box<str>),
    /// `.NAME.` (booleans are `.T.` and `.F.`).
    Enumeration(Name),
    /// `#123`.
    Reference(usize),
    /// `$`: no value.
    Unset,
    /// `*`: derived from other attributes.
    Derived,
    /// A list, its values kept by the exchange ([`Exchange::list`]).
    List(Stretch),
    /// A value of a named type, `NAME(value)` (e.g. `LENGTH_MEASURE(1.)`), its parameters
    /// kept by the exchange.
    Typed(Box<(Name, Stretch)>),
}

/// Where a list's values are among all the exchange's.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Stretch {
    pub(super) start: u32,
    pub(super) len: u32,
}

impl Stretch {
    pub(super) fn range(self) -> std::ops::Range<usize> {
        self.start as usize..(self.start + self.len) as usize
    }
}

impl Value {
    pub fn as_reference(&self) -> Option<usize> {
        match self {
            Value::Reference(r) => Some(*r),
            _ => None,
        }
    }

    pub fn as_integer(&self) -> Option<i64> {
        match self {
            Value::Integer(i) => Some(*i),
            Value::Real(x) if x.fract() == 0.0 => Some(*x as i64),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            Value::String(s) => Some(s),
            _ => None,
        }
    }

    /// `.T.` / `.F.` (`.U.`, unknown, counts as false).
    pub fn as_bool(&self) -> Option<bool> {
        match self {
            Value::Enumeration(e) => Some(&**e == "T"),
            _ => None,
        }
    }

    pub fn as_enumeration(&self) -> Option<&str> {
        match self {
            Value::Enumeration(e) => Some(e),
            _ => None,
        }
    }
}

/// A record's type name or an enumeration's value: shared by all its uses in a file
/// (a file has few distinct ones).
pub type Name = std::rc::Rc<str>;

/// An instance: one record `NAME(parameters)`, or several for a complex instance
/// `(A(...) B(...))` (its supertypes and subtypes, by name).
#[derive(Debug, Clone, Copy)]
pub struct Instance<'a> {
    x: &'a Exchange,
    records: &'a [(Name, Stretch)],
}

impl<'a> Instance<'a> {
    /// The type name of a simple instance (the first record's of a complex one).
    pub fn name(&self) -> &'a str {
        &self.records[0].0
    }

    pub fn is_complex(&self) -> bool {
        self.records.len() > 1
    }

    /// The records' type names and parameters.
    pub fn records(&self) -> impl Iterator<Item = (&'a str, &'a [Value])> + 'a {
        let x = self.x;
        self.records.iter().map(move |(n, v)| (&**n, &x.values[v.range()]))
    }

    /// The parameters of the record of type `name`.
    pub fn record(&self, name: &str) -> Option<&'a [Value]> {
        self.records.iter().find(|(n, _)| &**n == name).map(|(_, v)| &self.x.values[v.range()])
    }

    pub fn has(&self, name: &str) -> bool {
        self.records.iter().any(|(n, _)| &**n == name)
    }

    /// The parameters of a simple instance.
    pub fn params(&self) -> &'a [Value] {
        &self.x.values[self.records[0].1.range()]
    }
}

/// A file's header records and its data section's instances, by number.
#[derive(Debug, Default)]
pub struct Exchange {
    /// The header section's records (`FILE_DESCRIPTION`, `FILE_NAME`, `FILE_SCHEMA`, ...),
    /// in order.
    pub(super) header: Vec<(Name, Stretch)>,
    /// The instances' numbers and records, in file order.
    pub(super) instances: Vec<(usize, Stretch)>,
    /// All instances' records: type names and parameters.
    pub(super) records: Vec<(Name, Stretch)>,
    /// All lists' values.
    pub(super) values: Vec<Value>,
    /// Places in `instances` by number.
    pub(super) places: HashMap<usize, u32>,
    /// Places in `instances` by record type, in file order.
    pub(super) types: HashMap<Name, Vec<u32>>,
}

impl Exchange {
    pub fn get(&self, id: usize) -> Result<Instance<'_>> {
        match self.places.get(&id) {
            Some(&k) => Ok(self.instance(k)),
            None => Err(Error::InvalidInput(format!("STEP: no instance #{id}"))),
        }
    }

    fn instance(&self, place: u32) -> Instance<'_> {
        Instance { x: self, records: &self.records[self.instances[place as usize].1.range()] }
    }

    /// The instances' numbers in file order.
    pub fn ids(&self) -> impl Iterator<Item = usize> + '_ {
        self.instances.iter().map(|(id, _)| *id)
    }

    /// The instances of a type (simple ones, or complex ones with that record).
    pub fn of_type<'a>(&'a self, name: &str) -> impl Iterator<Item = (usize, Instance<'a>)> + 'a {
        let places = self.types.get(name).map_or(&[][..], Vec::as_slice);
        places.iter().map(|&k| (self.instances[k as usize].0, self.instance(k)))
    }

    /// The values of a list.
    pub fn list(&self, v: &Value) -> Option<&[Value]> {
        match v {
            Value::List(l) => Some(&self.values[l.range()]),
            _ => None,
        }
    }

    /// The type name and parameters of a typed value.
    pub fn typed<'a>(&'a self, v: &'a Value) -> Option<(&'a str, &'a [Value])> {
        match v {
            Value::Typed(t) => Some((&*t.0, &self.values[t.1.range()])),
            _ => None,
        }
    }

    /// A number, or a typed one (a measure).
    pub fn real(&self, v: &Value) -> Option<f64> {
        match v {
            Value::Real(x) => Some(*x),
            Value::Integer(i) => Some(*i as f64),
            Value::Typed(t) if t.1.len == 1 => self.real(&self.values[t.1.start as usize]),
            _ => None,
        }
    }
}

impl Exchange {
    /// The exchange as a STEP file's text: its header records and instances as read (or
    /// decoded), each on a line; reals in their shortest exact form.
    pub fn to_text(&self) -> String {
        let mut out = String::with_capacity(self.values.len() * 12);
        out.push_str("ISO-10303-21;\nHEADER;\n");
        for (name, params) in &self.header {
            self.write_record(&mut out, name, *params);
            out.push_str(";\n");
        }
        out.push_str("ENDSEC;\nDATA;\n");
        for &(id, records) in &self.instances {
            out.push('#');
            out.push_str(&id.to_string());
            out.push('=');
            let records = &self.records[records.range()];
            if records.len() == 1 {
                self.write_record(&mut out, &records[0].0, records[0].1);
            } else {
                out.push('(');
                for (name, params) in records {
                    self.write_record(&mut out, name, *params);
                }
                out.push(')');
            }
            out.push_str(";\n");
        }
        out.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
        out
    }

    fn write_record(&self, out: &mut String, name: &str, params: Stretch) {
        out.push_str(name);
        self.write_list(out, params);
    }

    fn write_list(&self, out: &mut String, list: Stretch) {
        out.push('(');
        for (k, v) in self.values[list.range()].iter().enumerate() {
            if k > 0 {
                out.push(',');
            }
            self.write_value(out, v);
        }
        out.push(')');
    }

    fn write_value(&self, out: &mut String, v: &Value) {
        use std::fmt::Write;
        match v {
            Value::Integer(i) => {
                let _ = write!(out, "{i}");
            }
            Value::Real(x) => out.push_str(&real(*x)),
            Value::String(s) => encode_string(out, s),
            Value::Binary(b) => {
                out.push('"');
                out.push_str(b);
                out.push('"');
            }
            Value::Enumeration(e) => {
                out.push('.');
                out.push_str(e);
                out.push('.');
            }
            Value::Reference(r) => {
                let _ = write!(out, "#{r}");
            }
            Value::Unset => out.push('$'),
            Value::Derived => out.push('*'),
            Value::List(l) => self.write_list(out, *l),
            Value::Typed(t) => self.write_record(out, &t.0, t.1),
        }
    }
}

/// A real in STEP's form: the shortest that reads back the same, with a point and `E` for
/// an exponent (`1.`, `1.5E-07`).
pub(super) fn real(v: f64) -> String {
    if !v.is_finite() {
        return "0.".into();
    }
    let s = format!("{v:?}");
    let (mantissa, exponent) = match s.split_once('e') {
        Some((m, e)) => (m.to_string(), Some(e.to_string())),
        None => (s, None),
    };
    let mantissa = if mantissa.contains('.') { mantissa } else { format!("{mantissa}.") };
    match exponent {
        Some(e) => format!("{mantissa}E{e}"),
        None => mantissa,
    }
}

/// A string in STEP's form: quoted, quotes and backslashes doubled, characters outside
/// printable ASCII as `\X2\....\X0\` (UTF-16 in hexadecimal), as [`decode`] reads them.
fn encode_string(out: &mut String, s: &str) {
    use std::fmt::Write;
    out.push('\'');
    let mut chars = s.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '\'' => out.push_str("''"),
            '\\' => out.push_str("\\\\"),
            ' '..='~' => out.push(c),
            _ => {
                // A run of them in one \X2\ ... \X0\.
                out.push_str("\\X2\\");
                let mut units = [0u16; 2];
                for u in c.encode_utf16(&mut units) {
                    let _ = write!(out, "{:04X}", u);
                }
                while let Some(&n) = chars.peek() {
                    if (' '..='~').contains(&n) {
                        break;
                    }
                    chars.next();
                    for u in n.encode_utf16(&mut units) {
                        let _ = write!(out, "{:04X}", u);
                    }
                }
                out.push_str("\\X0\\");
            }
        }
    }
    out.push('\'');
}

/// Parses a STEP file's text.
pub fn parse(text: &[u8]) -> Result<Exchange> {
    let mut p = Parser { s: text, i: 0, names: Default::default(), by_name: Vec::new(), stack: Vec::new(), values: Vec::new() };
    // Up to the data section.
    if text.len() >= u32::MAX as usize {
        return Err(Error::InvalidInput("STEP: files of 4 GB and more are not read".into()));
    }
    let mut exchange = Exchange::default();
    // The header's records, then the data section after it.
    let mut after_header = 0;
    if let Some(at) = find(text, b"HEADER;") {
        p.i = at + 7;
        loop {
            p.skip_space();
            if p.i >= p.s.len() {
                return Err(p.error("the end of the header"));
            }
            if p.s[p.i..].starts_with(b"ENDSEC") {
                p.i += 6;
                p.skip_space();
                p.expect(b';')?;
                break;
            }
            let (name, _) = p.keyword()?;
            p.skip_space();
            let params = p.list()?;
            p.skip_space();
            p.expect(b';')?;
            exchange.header.push((name, params));
        }
        after_header = p.i;
    }
    let data = find(&text[after_header..], b"DATA;").ok_or_else(|| Error::InvalidInput("STEP: no DATA section".into()))?;
    p.i = after_header + data + 5;
    // (About one instance per 60 bytes, and a value per 8, in files written by common
    // systems.)
    let rest = text.len() - p.i;
    exchange.instances.reserve(rest / 60);
    exchange.records.reserve(rest / 60);
    exchange.places.reserve(rest / 60);
    p.values.reserve(rest / 8);
    loop {
        p.skip_space();
        if p.i >= p.s.len() {
            break;
        }
        if p.s[p.i..].starts_with(b"ENDSEC") {
            break;
        }
        if p.s[p.i] != b'#' {
            return Err(p.error("an instance"));
        }
        p.i += 1;
        let id = p.integer()? as usize;
        p.skip_space();
        p.expect(b'=')?;
        p.skip_space();
        let place = exchange.instances.len() as u32;
        let first = exchange.records.len() as u32;
        let mut record = |p: &mut Parser| -> Result<()> {
            let (name, k) = p.keyword()?;
            let list = &mut p.by_name[k as usize];
            if list.last() != Some(&place) {
                list.push(place);
            }
            p.skip_space();
            let params = p.list()?;
            exchange.records.push((name, params));
            Ok(())
        };
        if p.peek() == Some(b'(') {
            // A complex instance: records in parentheses.
            p.i += 1;
            loop {
                p.skip_space();
                if p.peek() == Some(b')') {
                    p.i += 1;
                    break;
                }
                record(&mut p)?;
            }
        } else {
            record(&mut p)?;
        }
        p.skip_space();
        p.expect(b';')?;
        // (A number given twice: the last instance.)
        exchange.places.insert(id, place);
        let records = Stretch { start: first, len: exchange.records.len() as u32 - first };
        exchange.instances.push((id, records));
    }
    exchange.values = std::mem::take(&mut p.values);
    for (name, k) in p.names.into_values() {
        let places = std::mem::take(&mut p.by_name[k as usize]);
        if !places.is_empty() {
            exchange.types.insert(name, places);
        }
    }
    Ok(exchange)
}

fn find(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    haystack.windows(needle.len()).position(|w| w == needle)
}

struct Parser<'a> {
    s: &'a [u8],
    i: usize,
    /// The names met so far, each kept once, with their numbers.
    names: crate::hashing::FastMap<Box<[u8]>, (Name, u32)>,
    /// Places of the instances with records of each name, by the name's number.
    by_name: Vec<Vec<u32>>,
    /// The values of the lists being parsed, innermost last (each list's moved to `values`
    /// when it ends).
    stack: Vec<Value>,
    /// The values of the lists parsed.
    values: Vec<Value>,
}

impl Parser<'_> {
    fn error(&self, what: &str) -> Error {
        let line = self.s[..self.i.min(self.s.len())].iter().filter(|&&c| c == b'\n').count() + 1;
        Error::InvalidInput(format!("STEP: expected {what} at line {line}"))
    }

    fn peek(&self) -> Option<u8> {
        self.s.get(self.i).copied()
    }

    fn expect(&mut self, c: u8) -> Result<()> {
        if self.peek() == Some(c) {
            self.i += 1;
            Ok(())
        } else {
            Err(self.error(&format!("'{}'", c as char)))
        }
    }

    /// Whitespace and comments.
    fn skip_space(&mut self) {
        // Mostly there is none.
        if self.peek().is_some_and(|c| !c.is_ascii_whitespace() && c != b'/') {
            return;
        }
        loop {
            while self.i < self.s.len() && self.s[self.i].is_ascii_whitespace() {
                self.i += 1;
            }
            if self.s[self.i..].starts_with(b"/*") {
                match find(&self.s[self.i + 2..], b"*/") {
                    Some(end) => self.i += end + 4,
                    None => self.i = self.s.len(),
                }
            } else {
                break;
            }
        }
    }

    fn integer(&mut self) -> Result<i64> {
        let start = self.i;
        if matches!(self.peek(), Some(b'-') | Some(b'+')) {
            self.i += 1;
        }
        let digits = self.i;
        let mut value: i64 = 0;
        while let Some(c) = self.peek().filter(u8::is_ascii_digit) {
            value = value.checked_mul(10).and_then(|v| v.checked_add((c - b'0') as i64)).ok_or_else(|| self.error("an integer"))?;
            self.i += 1;
        }
        if digits == self.i {
            return Err(self.error("an integer"));
        }
        Ok(if self.s[start] == b'-' { -value } else { value })
    }

    /// A keyword, with its name's number.
    fn keyword(&mut self) -> Result<(Name, u32)> {
        let start = self.i;
        while matches!(self.peek(), Some(c) if c.is_ascii_alphanumeric() || c == b'_' || c == b'!') {
            self.i += 1;
        }
        if start == self.i {
            return Err(self.error("a keyword"));
        }
        Ok(self.name(start, self.i))
    }

    /// The name spelled by `s[start..end]`, shared with its other uses, and its number.
    fn name(&mut self, start: usize, end: usize) -> (Name, u32) {
        let bytes = &self.s[start..end];
        if let Some((name, k)) = self.names.get(bytes) {
            return (name.clone(), *k);
        }
        let name: Name = String::from_utf8_lossy(bytes).into();
        let k = self.by_name.len() as u32;
        self.by_name.push(Vec::new());
        self.names.insert(bytes.into(), (name.clone(), k));
        (name, k)
    }

    /// `( value, ... )`.
    fn list(&mut self) -> Result<Stretch> {
        self.expect(b'(')?;
        let start = self.stack.len();
        loop {
            self.skip_space();
            if self.peek() == Some(b')') {
                self.i += 1;
                let stretch = Stretch { start: self.values.len() as u32, len: (self.stack.len() - start) as u32 };
                self.values.extend(self.stack.drain(start..));
                return Ok(stretch);
            }
            let value = self.value()?;
            self.stack.push(value);
            self.skip_space();
            match self.peek() {
                Some(b',') => self.i += 1,
                Some(b')') => {}
                _ => return Err(self.error("',' or ')'")),
            }
        }
    }

    fn value(&mut self) -> Result<Value> {
        self.skip_space();
        match self.peek() {
            Some(b'#') => {
                self.i += 1;
                Ok(Value::Reference(self.integer()? as usize))
            }
            Some(b'$') => {
                self.i += 1;
                Ok(Value::Unset)
            }
            Some(b'*') => {
                self.i += 1;
                Ok(Value::Derived)
            }
            Some(b'(') => Ok(Value::List(self.list()?)),
            Some(b'\'') => self.string(),
            Some(b'.') => {
                self.i += 1;
                let start = self.i;
                while matches!(self.peek(), Some(c) if c != b'.') {
                    self.i += 1;
                }
                let (name, _) = self.name(start, self.i);
                self.expect(b'.')?;
                Ok(Value::Enumeration(name))
            }
            Some(b'"') => {
                // Binary: its digits.
                self.i += 1;
                let start = self.i;
                while matches!(self.peek(), Some(c) if c != b'"') {
                    self.i += 1;
                }
                let text = String::from_utf8_lossy(&self.s[start..self.i]).into();
                self.expect(b'"')?;
                Ok(Value::Binary(text))
            }
            Some(c) if c == b'-' || c == b'+' || c.is_ascii_digit() => self.number(),
            Some(c) if c.is_ascii_alphabetic() => {
                let (name, _) = self.keyword()?;
                self.skip_space();
                let params = self.list()?;
                Ok(Value::Typed(Box::new((name, params))))
            }
            _ => Err(self.error("a value")),
        }
    }

    fn number(&mut self) -> Result<Value> {
        let start = self.i;
        if matches!(self.peek(), Some(b'-') | Some(b'+')) {
            self.i += 1;
        }
        let mut real = false;
        while let Some(c) = self.peek() {
            if c.is_ascii_digit() {
                self.i += 1;
            } else if c == b'.' || c == b'E' || c == b'e' {
                real = true;
                self.i += 1;
                if (c == b'E' || c == b'e') && matches!(self.peek(), Some(b'-') | Some(b'+')) {
                    self.i += 1;
                }
            } else {
                break;
            }
        }
        let text = std::str::from_utf8(&self.s[start..self.i]).map_err(|_| self.error("a number"))?;
        if real {
            // ("1." and "1.E-07", valid STEP reals, parse as they are.)
            text.parse().map(Value::Real).map_err(|_| self.error("a number"))
        } else {
            text.parse().map(Value::Integer).map_err(|_| self.error("a number"))
        }
    }

    fn string(&mut self) -> Result<Value> {
        self.expect(b'\'')?;
        let mut out = Vec::new();
        loop {
            match self.peek() {
                None => return Err(self.error("the end of a string")),
                Some(b'\'') => {
                    if self.s.get(self.i + 1) == Some(&b'\'') {
                        out.push(b'\'');
                        self.i += 2;
                    } else {
                        self.i += 1;
                        break;
                    }
                }
                Some(c) => {
                    out.push(c);
                    self.i += 1;
                }
            }
        }
        Ok(Value::String(decode(&out).into()))
    }
}

/// A string's text: ISO 10303-21's `\X2\....\X0\` (UTF-16 hex) and `\X\..` (8-bit) escapes
/// decoded, other bytes as UTF-8.
fn decode(raw: &[u8]) -> String {
    let text = String::from_utf8_lossy(raw);
    if !text.contains('\\') {
        return text.into_owned();
    }
    let mut out = String::new();
    let mut rest: &str = &text;
    while let Some(k) = rest.find('\\') {
        out.push_str(&rest[..k]);
        rest = &rest[k..];
        if let Some(body) = rest.strip_prefix("\\X2\\") {
            let end = body.find("\\X0\\").unwrap_or(body.len());
            let units: Vec<u16> = (0..end / 4).filter_map(|i| u16::from_str_radix(&body[4 * i..4 * i + 4], 16).ok()).collect();
            out.push_str(&String::from_utf16_lossy(&units));
            rest = &body[(end + 4).min(body.len())..];
        } else if let Some(body) = rest.strip_prefix("\\X\\") {
            if let Some(b) = body.get(..2).and_then(|h| u8::from_str_radix(h, 16).ok()) {
                out.push(b as char);
                rest = &body[2..];
            } else {
                out.push('\\');
                rest = &rest[1..];
            }
        } else if let Some(body) = rest.strip_prefix("\\\\") {
            out.push('\\');
            rest = body;
        } else {
            out.push('\\');
            rest = &rest[1..];
        }
    }
    out.push_str(rest);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn values_are_small() {
        assert!(std::mem::size_of::<Value>() <= 24);
    }

    #[test]
    fn simple_and_complex_instances() {
        let text = b"ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n#1 = CARTESIAN_POINT('it''s',(0.,1.5,-2.E-03));\n#2 = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );\n#3 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#2,'distance_accuracy_value','');\n#4 = DIRECTION('',(0.,0.,1.)) /* a comment */;\nENDSEC;\nEND-ISO-10303-21;\n";
        let x = parse(text).unwrap();
        let p = x.get(1).unwrap();
        assert_eq!(p.name(), "CARTESIAN_POINT");
        assert_eq!(p.params()[0].as_str(), Some("it's"));
        let coords: Vec<f64> = x.list(&p.params()[1]).unwrap().iter().map(|v| x.real(v).unwrap()).collect();
        assert_eq!(coords, vec![0.0, 1.5, -2e-3]);
        let unit = x.get(2).unwrap();
        assert!(unit.is_complex());
        assert_eq!(unit.record("SI_UNIT").unwrap()[0].as_enumeration(), Some("MILLI"));
        assert_eq!(x.real(&x.get(3).unwrap().params()[0]), Some(1e-7));
        assert_eq!(x.ids().collect::<Vec<_>>(), vec![1, 2, 3, 4]);
        assert_eq!(x.of_type("NAMED_UNIT").map(|(id, _)| id).collect::<Vec<_>>(), vec![2]);
        assert_eq!(x.of_type("DIRECTION").map(|(id, _)| id).collect::<Vec<_>>(), vec![4]);
    }

    #[test]
    fn escaped_text() {
        assert_eq!(decode(b"caf\\X2\\00E9\\X0\\"), "caf\u{e9}");
    }
}
