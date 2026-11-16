"""OCCT's XCAF documents in its XML format (XmlXCAF), as CadQuery writes and reads them: an
assembly's labels (shapes, components placed, names, colours, layers, sub-shapes) and its
shapes, in OCCT's text BRep format.

Written as CadQuery's `toCAF` builds the document (label by label, in its order), read as its
`_importDoc` walks one."""

from __future__ import annotations

import datetime
import re
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape

from ..geom import Color, Location
from ..shapes import Compound, Shape
from ._occt_colors import NAMES as _COLOR_NAMES, SRGB as _COLOR_SRGB

# OCAF's GUIDs: tree nodes linking documents' parts, a label marked an assembly, names.
_DOC_TREE = "efd212eb-6dfd-11d4-b9c8-0060b0ee281b"
_SHAPE_REF = "5b896afe-3adf-11d4-b9b7-0060b0ee281b"
_ASSEMBLY = "5b896b00-3adf-11d4-b9b7-0060b0ee281b"
_COLOR_GEN = "efd212e4-6dfd-11d4-b9c8-0060b0ee281b"
_COLOR_SURF = "efd212e5-6dfd-11d4-b9c8-0060b0ee281b"
_LAYER_REF = "efd212e8-6dfd-11d4-b9c8-0060b0ee281b"
_NAME = "2a96b608-ec8b-11d0-bee7-080009dc3333"
_NS = "http://www.opencascade.org/OCAF/XML"
_XCAF_NS = "http://www.opencascade.org/OCAF/XML/XCAF"


def _nearest_color(rgba) -> int:
    """The index of OCCT's named colour nearest these sRGB components (as
    `Quantity_Color::Name`: squared distance in sRGB, single precision, the first nearest)."""
    import struct

    f = lambda x: struct.unpack("f", struct.pack("f", x))[0]
    r, g, b = (f(c) for c in rgba[:3])
    best, distance = 0, float("inf")
    for i, (x, y, z) in enumerate(_COLOR_SRGB):
        d = f(f(f(r - x) ** 2) + f(f(g - y) ** 2) + f(f(b - z) ** 2))
        if d < distance:
            best, distance = i, d
    return best


class _Node:
    """A tree node (or graph node) attribute: its label, father and children."""

    def __init__(self, label, tree: str, graph: bool = False):
        self.label, self.tree, self.graph = label, tree, graph
        self.fathers: list = []
        self.children: list = []
        self.id = 0


class _Label:
    def __init__(self, parent, tag: int):
        self.parent, self.tag = parent, tag
        self.children: dict = {}
        # Attributes in order: (kind, value).
        self.attributes: list = []
        self.next_tag = 0
        self.shape_ref = None

    def child(self, tag=None) -> "_Label":
        if tag is None:
            self.next_tag += 1
            tag = self.next_tag
        label = _Label(self, tag)
        self.children[tag] = label
        return label

    def get(self, kind):
        return next((v for k, v in self.attributes if k == kind), None)

    def set(self, kind, value) -> None:
        for i, (k, _) in enumerate(self.attributes):
            if k == kind:
                self.attributes[i] = (kind, value)
                return
        self.attributes.append((kind, value))

    def node(self, tree: str, graph: bool = False) -> _Node:
        """The label's tree (or graph) node of `tree`, made if it has none."""
        kind = ("graph" if graph else "tree", tree)
        found = self.get(kind)
        if found is None:
            found = _Node(self, tree, graph)
            self.attributes.append((kind, found))
        return found

    def entry(self) -> str:
        tags = []
        label = self
        while label is not None:
            tags.append(str(label.tag))
            label = label.parent
        return ":".join(reversed(tags))


def _link(child: _Node, father: _Node) -> None:
    child.fathers.append(father)
    father.children.append(child)


class _Document:
    """An XCAF document's labels: the shapes, colours and layers tools'."""

    def __init__(self):
        self.root = _Label(None, 0)
        main = self.root.child(1)
        main.attributes.append(("raw", "XCAFDoc_DocumentTool"))
        _link(main.node(_DOC_TREE), self.root.node(_DOC_TREE))
        self.tools = {}
        for tag, name, tool in [(1, "Shapes", "XCAFDoc_ShapeTool"), (2, "Colors", "XCAFDoc_ColorTool"), (3, "Layers", "XCAFDoc_LayerTool"), (4, "D&GTs", "XCAFDoc_DimTolTool"), (5, "Materials", "XCAFDoc_MaterialTool"), (7, "Views", "XCAFDoc_ViewTool"), (8, "Clipping Planes", "XCAFDoc_ClippingPlaneTool"), (9, "Notes", "XCAFDoc_NotesTool")]:
            label = main.child(tag)
            label.set("name", name)
            label.attributes.append(("raw", tool))
            self.tools[name] = label
        self.shapes = self.tools["Shapes"]
        self.colors: list = []
        self.layers: dict = {}

    # The shape tool's operations.

    def new_shape(self) -> _Label:
        return self.shapes.child()

    def add_component(self, assembly: _Label, referred: _Label, loc) -> _Label:
        if assembly.get("uattr") is None:
            assembly.set("uattr", _ASSEMBLY)
        component = assembly.child()
        component.set("location", loc)
        _link(component.node(_SHAPE_REF), referred.node(_SHAPE_REF))
        return component

    def add_subshape(self, label: _Label, shape) -> _Label:
        sub = label.child()
        sub.set("shape", shape)
        return sub

    def set_color(self, label: _Label, color, tree: str) -> None:
        rgba = tuple(color.toTuple())
        found = next((c for c, value in self.colors if all(abs(a - b) <= 1e-4 for a, b in zip(value, rgba))), None)
        if found is None:
            found = self.tools["Colors"].child()
            index = _nearest_color(rgba)
            hexa = "".join(f"{min(255, max(0, int(c * 255.0 + 0.5))):02X}" for c in rgba)
            found.set("color", index)
            found.set("name", f"{_COLOR_NAMES[index]} (#{hexa})")
            self.colors.append((found, rgba))
        _link(label.node(tree), found.node(tree))

    def set_layer(self, label: _Label, name: str) -> None:
        layer = self.layers.get(name)
        if layer is None:
            layer = self.tools["Layers"].child()
            layer.set("name", name)
            self.layers[name] = layer
        _link(label.node(_LAYER_REF, graph=True), layer.node(_LAYER_REF, graph=True))


def _is_identity(loc) -> bool:
    return all(abs(v - (1.0 if i == j else 0.0)) <= 1e-12 for i, row in enumerate(loc.rows) for j, v in enumerate(row))


def _build(assy) -> _Document:
    """The document CadQuery's `toCAF` makes of the assembly (colours not as STEP's)."""
    doc = _Document()
    unique: dict = {}
    compounds: dict = {}

    def visit(el, ancestor):
        subassy = None
        if el.children:
            subassy = doc.new_shape()
            subassy.set("name", el.name)
        color = el.color if el.color else None
        lab = None
        if el.obj:
            key0 = (color.toTuple() if color is not None else None, id(el.obj))
            if key0 in unique:
                lab = unique[key0]
            else:
                lab = doc.new_shape()
                compound = compounds.get(id(el.obj))
                if compound is None:
                    compound = Compound.makeCompound(el.shapes)
                    compounds[id(el.obj)] = compound
                lab.set("shape", compound)
                lab.set("name", f"{el.name}_part" if el.children else el.name)
                unique[key0] = lab
                if getattr(el, "material", None):
                    raise NotImplementedError("materials in OCAF XML documents")
            marks: dict = {}
            for attribute, slot in (("_subshape_names", 0), ("_subshape_colors", 1), ("_subshape_layers", 2)):
                for s, value in getattr(el, attribute, {}).items():
                    marks.setdefault(s, [None, None, None])[slot] = value
            for s, (name, sub_color, layer) in marks.items():
                sub = doc.add_subshape(lab, s)
                if name is not None:
                    sub.set("name", name)
                if sub_color is not None:
                    doc.set_color(sub, sub_color, _COLOR_GEN)
                if layer is not None:
                    doc.set_layer(sub, layer)
            if el.children:
                lab = doc.add_component(subassy, lab, Location())
                lab.set("name", f"{el.name}_part")
            elif ancestor is not None:
                lab = doc.add_component(ancestor, lab, el.loc)
                lab.set("name", el.name)
        if color is not None:
            if el.children:
                doc.set_color(subassy, color, _COLOR_SURF)
            if el.obj:
                doc.set_color(lab, color, _COLOR_SURF)
        for child in el.children:
            visit(child, subassy)
        if ancestor is not None and el.children:
            doc.add_component(ancestor, subassy, el.loc)
        elif ancestor is not None:
            pass
        elif el.children:
            if not _is_identity(assy.loc):
                # (A new free label referring to the top, placed: the shape tool's SetLocation.)
                top = doc.new_shape()
                top.set("location", assy.loc)
                _link(top.node(_SHAPE_REF), subassy.node(_SHAPE_REF))
                top.set("name", assy.name)
        elif el.obj:
            rv = doc.new_shape()
            component = doc.add_component(rv, lab, el.loc)
            component.set("name", el.name)
        else:
            raise ValueError("Cannot convert an empty assembly to CAF")

    visit(assy, None)
    # Each assembly's shape: its components' shapes placed, in a compound (deepest first).
    def shape_of(label: _Label):
        if label.get("uattr") is None:
            return label.get("shape")
        if label.get("shape") is None:
            parts = []
            for component in label.children.values():
                ref = component.get(("tree", _SHAPE_REF))
                if ref is None or not ref.fathers:
                    continue
                loc = component.get("location")
                s = shape_of(ref.fathers[0].label)
                parts.append(s if loc is None or _is_identity(loc) else s.moved(loc))
            label.set("shape", Compound.makeCompound(parts))
            # (The shape first among the label's attributes, as OCCT writes it.)
            label.attributes.sort(key=lambda a: a[0] != "shape")
        return label.get("shape")

    for label in list(doc.shapes.children.values()):
        shape_of(label)
    return doc


def _locations(lines: list) -> list:
    """The Locations table of a text BRep's lines: each entry's 3 x 4 rows (index 1 first)."""
    i = next(k for k, line in enumerate(lines) if line.startswith("Locations"))
    count = int(lines[i].split()[1])
    out: list = []
    k = i + 1
    identity = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]

    def compose(a, b):
        return [[sum(a[r][m] * (b[m][c] if m < 3 else 0.0) for m in range(3)) + (a[r][3] if c == 3 else 0.0) for c in range(4)] for r in range(3)]

    def invert(a):
        rt = [[a[c][r] for c in range(3)] for r in range(3)]
        t = [-sum(rt[r][m] * a[m][3] for m in range(3)) for r in range(3)]
        return [rt[r] + [t[r]] for r in range(3)]

    for _ in range(count):
        kind = lines[k].split()[0]
        if kind == "1":
            rows = [[float(v) for v in lines[k + 1 + r].split()] for r in range(3)]
            out.append(rows)
            k += 4
        else:
            numbers = []
            while True:
                numbers += [int(v) for v in lines[k].split()[(1 if not numbers else 0):]]
                k += 1
                if numbers and numbers[-1] == 0 and len(numbers) % 2 == 1:
                    break
            # (Each item's power multiplied in on the left of those before, as OCCT reads.)
            product = identity
            for index, power in zip(numbers[0:-1:2], numbers[1:-1:2]):
                m = out[index - 1]
                for _ in range(abs(power)):
                    product = compose(m if power > 0 else invert(m), product)
            out.append(product)
    return out


def _shape_text(doc: _Document):
    """The document's shapes in one text BRep (each label's shape and sub-shape a child of a
    compound written last), and each such shape's (sign, index, location) there."""
    labels = []

    def walk(label):
        if label.get("shape") is not None:
            labels.append(label)
        for child in label.children.values():
            walk(child)

    walk(doc.shapes)
    # (And, for each placement a label holds, the shape it refers to so placed: so that the
    # placement is in the locations written, as the top level's own.)
    extra = []

    def placed(label):
        loc = label.get("location")
        ref = label.get(("tree", _SHAPE_REF))
        if loc is not None and not _is_identity(loc) and ref is not None and ref.fathers:
            extra.append(ref.fathers[0].label.get("shape").moved(loc))
        for child in label.children.values():
            placed(child)

    placed(doc.shapes)
    whole = Compound.makeCompound([label.get("shape") for label in labels] + extra)
    text = whole.__getstate__()[0].decode()
    lines = text.splitlines()
    t = next(k for k, line in enumerate(lines) if line.startswith("TShapes"))
    count = int(lines[t].split()[1])
    # The compound's sub-shapes: the last entry's, before the root line.
    end = max(k for k, line in enumerate(lines) if line.strip())
    refs = []
    k = end - 1
    while not lines[k].strip().endswith("*"):
        k -= 1
    tokens = lines[k].split()
    for ref, loc in zip(tokens[0:-1:2], tokens[1:-1:2]):
        refs.append((ref[0], count - int(ref[1:]) + 1, int(loc)))
    if len(refs) != len(labels) + len(extra):
        raise RuntimeError("the shapes of the document were not written one each")
    for label, ref in zip(labels, refs):
        label.shape_ref = ref
    body = "\n".join(lines[:end]) + "\n"
    return body, _locations(lines)


def write_xml(assy, path: str) -> bool:
    doc = _build(assy)
    body, table = _shape_text(doc)

    def location_id(loc) -> int:
        if loc is None or _is_identity(loc):
            return 0
        for i, rows in enumerate(table):
            if all(abs(a - b) <= 1e-9 * (1.0 + abs(b)) for r, q in zip(rows, loc.rows) for a, b in zip(r, q)):
                return i + 1
        raise RuntimeError("a component's placement missing from the document's shapes")

    ids = [0]
    out = []
    counter = [0]

    def write(label, depth):
        pad = " " * depth
        out.append(f'{pad}<label tag="{label.tag}">')
        inner = " " * (depth + 1)
        for kind, value in label.attributes:
            counter[0] += 1
            n = counter[0]
            key = kind[0] if isinstance(kind, tuple) else kind
            if key == "raw":
                out.append(f'{inner}<{value} id="{n}"/>')
            elif key == "name":
                out.append(f'{inner}<TDataStd_Name id="{n}" nameguid="{_NAME}">{escape(value)}</TDataStd_Name>')
            elif key == "uattr":
                out.append(f'{inner}<TDataStd_UAttribute id="{n}" guid="{value}"/>')
            elif key == "color":
                out.append(f'{inner}<xcaf:Color id="{n}">{value}</xcaf:Color>')
            elif key == "location":
                loc = location_id(value)
                if loc:
                    out.append(f'{inner}<xcaf:Location id="{n}">\n{inner} <location locId="{loc}"/>\n{inner}</xcaf:Location>')
                else:
                    out.append(f'{inner}<xcaf:Location id="{n}"/>')
            elif key == "shape":
                sign, index, loc = label.shape_ref
                top = label.parent is doc.shapes
                located = f' location="{loc}"' if loc else ""
                version = ' version="1"' if top else ""
                out.append(f'{inner}<TNaming_NamedShape id="{n}" evolution="primitive"{version}>\n{inner} <olds/>\n{inner} <news>\n{inner}  <shape tshape="{sign}{index}"{located} index="1"/>\n{inner} </news>\n{inner}</TNaming_NamedShape>')
            elif key == "tree":
                children = f' children="{"".join(f"{c.id} " for c in value.children)}"' if value.children else ""
                out.append(f'{inner}<TDataStd_TreeNode id="{n}" treeid="{value.tree}"{children}/>')
            elif key == "graph":
                children = f' children="{"".join(f"{c.id} " for c in value.children)}"' if value.children else ""
                fathers = f' fathers="{"".join(f"{c.id} " for c in value.fathers)}"' if value.fathers else ""
                out.append(f'{inner}<xcaf:GraphNode id="{n}" treeid="{value.tree}"{children}{fathers}/>')
            if key not in ("raw", "name", "uattr", "color", "location", "shape", "tree", "graph"):
                raise RuntimeError(f"unknown attribute {kind}")
        if label.next_tag:
            counter[0] += 1
            out.append(f'{inner}<TDF_TagSource id="{counter[0]}">{label.next_tag}</TDF_TagSource>')
        for child in label.children.values():
            write(child, depth + 1)
        out.append(f"{pad}</label>")

    def renumber(label):
        # Ids in writing order (tree nodes' children are named by them): a label's
        # attributes, its tag source, then its children's.
        for kind, value in label.attributes:
            ids[0] += 1
            if isinstance(value, _Node):
                value.id = ids[0]
        if label.next_tag:
            ids[0] += 1
        for child in label.children.values():
            renumber(child)

    renumber(doc.root)
    total = ids[0]
    write(doc.root, 1)
    header = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<document format="XmlXCAF" xmlns="{_NS}" xmlns:xcaf="{_XCAF_NS}" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="{_NS}">\n'
        f' <info date="{datetime.date.today().isoformat()}" schemav="0" DocVersion="12" objnb="{total}">\n'
        "  <iitem>Copyright: Open Cascade, 2001-2002</iitem>\n"
        "  <iitem>REFERENCE_COUNTER: 0</iitem>\n"
        "  <iitem>MODIFICATION_COUNTER: 1</iitem>\n"
        " </info>\n"
        " <comments/>\n"
    )
    with open(path, "w", encoding="utf-8") as f:
        f.write(header + "\n".join(out) + "\n <shapes>\n" + body + "\n </shapes>\n</document>\n")
    return True


# Reading.


class _Read:
    """A document read: its labels (by entry), attributes and shapes."""

    def __init__(self, path: str):
        tree = ET.parse(path)
        root = tree.getroot()
        if not root.tag.endswith("document"):
            raise ValueError(f"not an OCAF XML document: {path}")
        self.labels: dict = {}
        self.nodes: dict = {}
        self.text = ""
        for child in root:
            name = child.tag.split("}")[-1]
            if name == "label":
                self._label(child, None)
            elif name == "shapes":
                self.text = child.text or ""
        # Tree nodes' fathers, from their fathers' children; graph nodes' both ways.
        for node in self.nodes.values():
            for c in node["children"]:
                if c in self.nodes:
                    self.nodes[c]["father"] = node
        self.shapes_cache: dict = {}
        self._shape_lines = self.text.splitlines()
        self._table = _locations(self._shape_lines) if any(l.startswith("Locations") for l in self._shape_lines) else []

    def _label(self, element, parent):
        tag = int(element.get("tag"))
        entry = f"{parent}:{tag}" if parent is not None else str(tag)
        label = {"entry": entry, "tag": tag, "parent": parent, "children": [], "name": None, "shape": None, "location": None, "uattr": set(), "trees": {}, "graphs": {}, "color": None}
        self.labels[entry] = label
        if parent is not None:
            self.labels[parent]["children"].append(entry)
        for a in element:
            kind = a.tag.split("}")[-1]
            if kind == "label":
                self._label(a, entry)
            elif kind == "TDataStd_Name":
                label["name"] = a.text or ""
            elif kind == "TNaming_NamedShape":
                s = a.find(f"{{{_NS}}}news/{{{_NS}}}shape")
                if s is not None:
                    t = s.get("tshape")
                    label["shape"] = (t[0], int(t[1:]), int(s.get("location") or 0))
            elif kind == "Location":
                loc = a.find(f"{{{_NS}}}location")
                label["location"] = int(loc.get("locId")) if loc is not None else 0
            elif kind == "TDataStd_UAttribute":
                label["uattr"].add(a.get("guid"))
            elif kind in ("TDataStd_TreeNode", "GraphNode"):
                node = {"id": a.get("id"), "label": entry, "tree": a.get("treeid"), "children": (a.get("children") or "").split(), "fathers": (a.get("fathers") or "").split(), "father": None}
                self.nodes[node["id"]] = node
                (label["graphs"] if kind == "GraphNode" else label["trees"])[node["tree"]] = node
            elif kind == "Color":
                label["color"] = (a.text or "").strip()

    # The shape tool's questions.

    def father(self, entry: str, tree: str):
        node = self.labels[entry]["trees"].get(tree)
        if node is None or node["father"] is None:
            return None
        return node["father"]["label"]

    def is_assembly(self, entry: str) -> bool:
        return _ASSEMBLY in self.labels[entry]["uattr"]

    def is_reference(self, entry: str) -> bool:
        return self.father(entry, _SHAPE_REF) is not None

    def components(self, entry: str) -> list:
        if not self.is_assembly(entry):
            return []
        return [c for c in self.labels[entry]["children"] if self.is_reference(c)]

    def location(self, entry: str) -> Location:
        loc = self.labels[entry]["location"]
        if not loc:
            return Location()
        return Location([list(r) for r in self._table[loc - 1]])

    def shapes_root(self) -> str:
        return next(e for e, l in self.labels.items() if l["name"] == "Shapes" and l["parent"] is not None and self.labels[l["parent"]]["parent"] == "0")

    def free_shapes(self) -> list:
        root = self.shapes_root()
        referred = {self.father(e, _SHAPE_REF) for e in self.labels if self.is_reference(e)}
        return [c for c in self.labels[root]["children"] if c not in referred and (self.labels[c]["shape"] is not None or self.is_reference(c) or self.is_assembly(c))]

    def color(self, entry: str, trees=(_COLOR_SURF, _COLOR_GEN)):
        for tree in trees:
            father = self.father(entry, tree)
            if father is not None:
                value = self.labels[father]["color"]
                return _color_of(value)
        return None

    def layers(self, entry: str) -> list:
        node = self.labels[entry]["graphs"].get(_LAYER_REF)
        if node is None:
            return []
        return [self.labels[self.nodes[f]["label"]]["name"] for f in node["fathers"] if f in self.nodes]

    def shapes(self, entries: list) -> list:
        """The shapes of these labels: the text BRep read once with a compound of them all
        made its root."""
        wanted = [self.labels[e]["shape"] for e in entries if self.labels[e]["shape"] is not None and e not in self.shapes_cache]
        if wanted:
            lines = list(self._shape_lines)
            t = next(k for k, line in enumerate(lines) if line.startswith("TShapes"))
            count = int(lines[t].split()[1])
            lines[t] = f"TShapes {count + 1}"
            # (Sub-shapes are named counting back from the last shape: each one further back
            # by the shape added after them. Each shape's sub-shapes end it, on a line ending
            # with "*".)
            for k in range(t + 1, len(lines)):
                tokens = lines[k].split()
                if tokens and tokens[-1] == "*" and len(tokens) % 2 == 1:
                    pairs = [f"{r[0]}{int(r[1:]) + 1} {l}" for r, l in zip(tokens[0:-1:2], tokens[1:-1:2])]
                    lines[k] = " ".join(pairs + ["*"])
            refs = " ".join(f"{sign}{count + 1 - index + 1} {loc}" for sign, index, loc in wanted)
            while lines and not lines[-1].strip():
                lines.pop()
            text = "\n".join(lines) + f"\nCo\n\n1100000\n{refs} *\n\n+1 0\n"
            s = Shape.__new__(Shape)
            s.__setstate__((text.encode(), False))
            whole = Shape.cast(s.wrapped)
            pieces = list(whole)
            todo = [e for e in entries if self.labels[e]["shape"] is not None and e not in self.shapes_cache]
            for e, piece in zip(todo, pieces):
                self.shapes_cache[e] = piece
        return [self.shapes_cache.get(e) for e in entries]

    def shape(self, entry: str):
        return self.shapes([entry])[0]


def _within(s, whole):
    """The sub-shape `s` as `whole` holds it: a wire referred to on its own reads apart from
    its face's (a kentron wire is its face's, its edges used on it), so the face's wire with
    the same edges."""
    from .shapes import isSubshape

    if s.ShapeType() != "Wire" or isSubshape(s, whole):
        return s
    edges = s.Edges()
    for w in whole.Wires():
        mine = w.Edges()
        if len(mine) == len(edges) and all(any(e.isSame(m) for m in mine) for e in edges):
            return w
    return s


def _color_of(value: str):
    """A colour attribute's value: a named colour's index (OCCT's XML), else components."""
    parts = value.split()
    if len(parts) == 1 and re.fullmatch(r"-?\d+", parts[0]):
        r, g, b = _COLOR_SRGB[int(parts[0])]
        return Color(r, g, b, 1.0)
    numbers = [float(p) for p in parts]
    return Color(*numbers[:4], False) if len(numbers) >= 3 else None


def read_xml(assy, path: str) -> None:
    """Fills `assy` from an OCAF XML document, as CadQuery's `_importDoc` does."""
    doc = _Read(path)

    def process(entry: str, parent):
        for comp in doc.components(entry):
            comp_name = doc.labels[comp]["name"] or ""
            loc = doc.location(comp)
            ref = doc.father(comp, _SHAPE_REF)
            color = doc.color(comp)
            if doc.is_assembly(ref):
                ref_name = doc.labels[ref]["name"] or ""
                sub = assy.__class__(name=ref_name)
                process(ref, sub)
                parent.add(sub, loc=loc, name=ref_name, color=color)
            elif doc.labels[ref]["shape"] is not None:
                ref_name = comp_name
                shape = doc.shape(ref)
                if color is None:
                    color = doc.color(ref, (_COLOR_SURF,))
                if ref_name.endswith("_part"):
                    parent.obj, parent.loc, parent.color = shape, loc, color
                    current = parent
                else:
                    parent.add(assy.__class__(shape, loc=loc, name=comp_name, color=color))
                    current = parent[comp_name]
                subs = [c for c in doc.labels[ref]["children"] if doc.labels[c]["shape"] is not None]
                for sub, s in zip(subs, doc.shapes(subs)):
                    s = _within(s, shape)
                    name = doc.labels[sub]["name"]
                    if name:
                        current.addSubshape(s, name=name)
                    for layer in doc.layers(sub):
                        current.addSubshape(s, layer=layer)
                    sub_color = doc.color(sub)
                    if sub_color:
                        current.addSubshape(s, color=sub_color)
        return parent

    # (All the shapes read at once: a sub-shape then the very face or edge of its part.)
    doc.shapes([e for e, label in doc.labels.items() if label["shape"] is not None])
    free = doc.free_shapes()
    if not free:
        raise ValueError(f"no shapes in the document {path}")
    top = free[0]
    color = doc.color(top)
    if doc.is_reference(top):
        top = doc.father(top, _SHAPE_REF)
    if doc.labels[top]["parent"] == doc.shapes_root() and doc.is_assembly(top):
        assy.objects.pop(assy.name)
        assy.name = doc.labels[top]["name"] or ""
        assy.objects[assy.name] = assy
        if color:
            assy.color = color
        loc = doc.location(top)
        assy.loc = loc
        imported = assy.__class__()
        process(top, imported)
        if assy.name in imported:
            imported = imported[assy.name]
            assy.loc = imported.loc
        for child in imported.children:
            assy.add(child, name=child.name, color=child.color, loc=child.loc)
    else:
        raise ValueError("Step file does not contain an assembly")
