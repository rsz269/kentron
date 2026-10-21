"""CadQuery's string selector syntax, parsed by hand (without pyparsing): `>Z`, `|(1, 0, 0)`,
`%CIRCLE`, `<<X[-2]`, `front`, combined with `and`, `or`, `exc` (`except`) and `not`, and
parentheses.

It reads what CadQuery's pyparsing grammar (`cadquery.selectors._makeGrammar` and
`_makeExpressionGrammar`) reads, the way pyparsing reads it:

* each token is matched at the next non-whitespace character, with no word boundaries
  (`>Xand<Y` is `>X and <Y`); a number takes no whitespace inside (`-1.5`, `2.`, no exponent);
* of the atoms' forms, the first that matches wins (a bare direction, `%type`, `>`/`<` with
  a direction and an optional `[index]`, `>>`/`<<` likewise, `|`/`#`/`+`/`-` with a
  direction, a named view); of several words where one begins another, the longer;
* `and` binds tighter than `or`, `or` tighter than `exc`, and `not` least: `not` takes
  everything after it, to the end or the closing parenthesis (`not A and B` is
  `not (A and B)`, `A and not B or C` is `A and not (B or C)`).

`parse` builds the result with the callables it is given, as the grammar's parse actions do.
"""

from __future__ import annotations

from typing import Callable, Iterable, List, Optional, Tuple

_WHITESPACE = " \t\n\r"
_DIRECTIONS = ("X", "Y", "Z", "XY", "XZ", "YZ")
_VIEWS = ("front", "back", "left", "right", "top", "bottom")


class SelectorSyntaxError(ValueError):
    """A string that is no selector: where reading it stopped."""

    def __init__(self, text: str, at: int):
        super().__init__(f"invalid selector {text!r} (at column {at + 1})")
        self.text = text
        self.at = at


class Fields(dict):
    """An atom's parts by name, as pyparsing's results give them: `name in fields`, and
    `fields.name` (empty where absent)."""

    def __getattr__(self, name):
        return self.get(name, "")


class Index(list):
    """An index's parts (`["-", "2"]`), as pyparsing's group gives them."""

    def asList(self) -> list:
        return list(self)


class VectorText:
    """A vector's coordinates as written."""

    def __init__(self, x: str, y: str, z: str):
        self.x, self.y, self.z = x, y, z


class _Parser:
    def __init__(self, text: str, types: Iterable[str], atom, and_, or_, exc, not_):
        self.text = text
        # (Longest first: of words where one begins another, the longer is taken.)
        self.types = sorted({t.upper() for t in types}, key=len, reverse=True)
        self.atom_of, self.and_, self.or_, self.exc, self.not_ = atom, and_, or_, exc, not_
        # How far reading got, for the error.
        self.reached = 0

    def skip(self, i: int) -> int:
        while i < len(self.text) and self.text[i] in _WHITESPACE:
            i += 1
        self.reached = max(self.reached, i)
        return i

    def word(self, i: int, options: Iterable[str], caseless: bool = False) -> Optional[Tuple[str, int]]:
        """The longest of `options` at `i` (after whitespace): it and where it ends."""
        i = self.skip(i)
        found = None
        for option in options:
            piece = self.text[i:i + len(option)]
            if (piece.upper() == option.upper() if caseless else piece == option) and (found is None or len(option) > len(found)):
                found = option
        return None if found is None else (found, i + len(found))

    def digits(self, i: int) -> int:
        """Where the digits from `i` end (`i` itself if none)."""
        while i < len(self.text) and self.text[i].isdigit():
            i += 1
        return i

    def number(self, i: int) -> Optional[Tuple[str, int]]:
        """An optionally signed integer with an optional point and fraction, nothing between."""
        i = self.skip(i)
        j = i + 1 if self.text[i:i + 1] in ("+", "-") else i
        k = self.digits(j)
        if k == j:
            return None
        if self.text[k:k + 1] == ".":
            k = self.digits(k + 1)
        self.reached = max(self.reached, k)
        return self.text[i:k], k

    def vector(self, i: int) -> Optional[Tuple[VectorText, int]]:
        """`(x, y, z)`, whitespace allowed between its parts."""
        found = self.word(i, ("(",))
        if found is None:
            return None
        coordinates: List[str] = []
        at = found[1]
        for separator in (",", ",", ")"):
            n = self.number(at)
            if n is None:
                return None
            coordinates.append(n[0])
            s = self.word(n[1], (separator,))
            if s is None:
                return None
            at = s[1]
        return VectorText(*coordinates), at

    def direction(self, i: int, fields: Fields) -> Optional[int]:
        """A direction at `i` into `fields` (`simple_dir` or `vector_dir`): where it ends."""
        simple = self.word(i, _DIRECTIONS)
        if simple is not None:
            fields["simple_dir"] = simple[0]
            return simple[1]
        vector = self.vector(i)
        if vector is not None:
            fields["vector_dir"] = vector[0]
            return vector[1]
        return None

    def index(self, i: int) -> Optional[Tuple[Index, int]]:
        """`[n]` or `[-n]`, whitespace allowed between its parts."""
        found = self.word(i, ("[",))
        if found is None:
            return None
        parts = Index()
        at = found[1]
        minus = self.word(at, ("-",))
        if minus is not None:
            parts.append("-")
            at = minus[1]
        start = self.skip(at)
        end = self.digits(start)
        if end == start:
            return None
        parts.append(self.text[start:end])
        close = self.word(end, ("]",))
        if close is None:
            return None
        return parts, close[1]

    def atom(self, i: int):
        """A simple selector at `i`, the first of its forms that matches: it and where it ends."""
        fields = Fields()
        end = self.direction(i, fields)
        if end is not None:
            fields["only_dir"] = True
            return self.atom_of(fields), end
        fields = Fields()
        found = self.word(i, ("%",))
        if found is not None:
            kind = self.word(found[1], self.types, caseless=True)
            if kind is not None:
                fields["type_op"], fields["cq_type"] = "%", kind[0]
                return self.atom_of(fields), kind[1]
        for key, operators in (("dir_op", (">", "<")), ("center_nth_op", (">>", "<<"))):
            fields = Fields()
            found = self.word(i, operators)
            if found is not None:
                end = self.direction(found[1], fields)
                if end is not None:
                    fields[key] = found[0]
                    index = self.index(end)
                    if index is not None:
                        fields["index"], end = index
                    return self.atom_of(fields), end
        fields = Fields()
        found = self.word(i, ("|", "#", "+", "-"))
        if found is not None:
            end = self.direction(found[1], fields)
            if end is not None:
                fields["other_op"] = found[0]
                return self.atom_of(fields), end
        found = self.word(i, _VIEWS)
        if found is not None:
            return self.atom_of(Fields(named_view=found[0])), found[1]
        return None

    def operand(self, i: int):
        """An atom, an expression in parentheses, or `not` before an expression (to as far
        as it reaches)."""
        found = self.atom(i)
        if found is not None:
            return found
        negated = self.word(i, ("not",))
        if negated is not None:
            inner = self.expression(negated[1])
            if inner is not None:
                return self.not_(inner[0]), inner[1]
        open_ = self.word(i, ("(",))
        if open_ is None:
            return None
        inner = self.expression(open_[1])
        if inner is None:
            return None
        close = self.word(inner[1], (")",))
        if close is None:
            return None
        return inner[0], close[1]

    def chain(self, i: int, operators: Tuple[str, ...], lower: Callable, combine: Callable):
        """`lower` once or more, joined by one of `operators` (an operator whose right side
        does not follow is not taken)."""
        first = lower(i)
        if first is None:
            return None
        items, at = [first[0]], first[1]
        while True:
            op = self.word(at, operators)
            if op is None:
                break
            right = lower(op[1])
            if right is None:
                break
            items.append(right[0])
            at = right[1]
        return (combine(items) if len(items) > 1 else items[0]), at

    def conjunction(self, i: int):
        return self.chain(i, ("and",), self.operand, self.and_)

    def disjunction(self, i: int):
        return self.chain(i, ("or",), self.conjunction, self.or_)

    def difference(self, i: int):
        return self.chain(i, ("exc", "except"), self.disjunction, self.exc)

    def expression(self, i: int):
        """`not` (any number of times) before a difference."""
        found = self.word(i, ("not",))
        if found is not None:
            inner = self.expression(found[1])
            if inner is not None:
                return self.not_(inner[0]), inner[1]
        return self.difference(i)


def parse(
    text: str,
    types: Iterable[str],
    atom: Callable[[Fields], object],
    and_: Callable[[list], object],
    or_: Callable[[list], object],
    exc: Callable[[list], object],
    not_: Callable[[object], object],
):
    """The selector `text` describes, built with `atom` (from an atom's `Fields`), `and_`,
    `or_` and `exc` (from the operands of a chain of one operator, in order) and `not_`;
    `types` are the names `%` takes (in any case). Raises SelectorSyntaxError."""
    parser = _Parser(text, types, atom, and_, or_, exc, not_)
    found = parser.expression(0)
    if found is not None and parser.skip(found[1]) == len(text):
        return found[0]
    raise SelectorSyntaxError(text, min(parser.reached, len(text)))
