"""Assemblies: named, placed and coloured shapes in a tree (CadQuery's `Assembly`)."""

from __future__ import annotations

import uuid
from typing import Dict, Iterator, List, Optional, Tuple

from .geom import Color, Location
from .shapes import Compound, Shape


class Assembly:
    """A node of an assembly tree: an optional shape (or workplane), its placement, name and
    colour, and child assemblies."""

    def __init__(self, obj=None, loc: Optional[Location] = None, name: Optional[str] = None, color: Optional[Color] = None, metadata=None, material=None):
        self.obj = obj
        self.loc = loc if loc else Location()
        self.name = name if name else str(uuid.uuid4())
        self.color = color if color else None
        self.metadata = metadata if metadata else {}
        self.material = material
        self.parent: Optional[Assembly] = None
        self.children: List[Assembly] = []
        self.objects: Dict[str, Assembly] = {self.name: self}

    def _copy(self) -> "Assembly":
        rv = self.__class__(self.obj, self.loc, self.name, self.color, self.metadata, self.material)
        for ch in self.children:
            child = ch._copy()
            child.parent = rv
            rv.children.append(child)
            rv.objects.update(child._flatten())
        return rv

    def _flatten(self, parents=()) -> Dict[str, "Assembly"]:
        path = "/".join(list(parents) + [self.name])
        rv = {path: self}
        for ch in self.children:
            rv.update(ch._flatten(list(parents) + [self.name]))
        return rv

    def add(self, arg, loc: Optional[Location] = None, name: Optional[str] = None, color: Optional[Color] = None, metadata=None, material=None) -> "Assembly":
        if isinstance(arg, Assembly):
            name = name or arg.name
            if name in self.objects:
                raise ValueError(f"Unique name is required. {name} is already in the assembly")
            sub = arg._copy()
            sub.loc = loc if loc else arg.loc
            sub.name = name
            sub.color = color if color else arg.color
            sub.metadata = metadata if metadata else arg.metadata
            sub.parent = self
            self.children.append(sub)
            self.objects.update(sub._flatten())
        else:
            self.add(self.__class__(arg, loc, name, color, metadata, material))
        return self

    @property
    def shapes(self) -> List[Shape]:
        from .workplane import Workplane

        if isinstance(self.obj, Shape):
            return [self.obj]
        if isinstance(self.obj, Workplane):
            return [o for o in self.obj.vals() if isinstance(o, Shape)]
        return []

    def traverse(self) -> Iterator[Tuple[str, "Assembly"]]:
        for ch in self.children:
            yield from ch.traverse()
        yield (self.name, self)

    def __iter__(self, loc: Optional[Location] = None, name: Optional[str] = None, color: Optional[Color] = None):
        name = f"{name}/{self.name}" if name else self.name
        loc = loc * self.loc if loc else self.loc
        color = self.color if self.color else color
        if self.obj is not None:
            shapes = self.shapes
            yield (shapes[0] if isinstance(self.obj, Shape) else Compound.makeCompound(shapes)), name, loc, color
        for ch in self.children:
            yield from ch.__iter__(loc, name, color)

    def _step_tree(self):
        """The assembly as nested `(name, shape, placement rows, rgb, children)` tuples."""
        shapes = self.shapes
        shape = None
        if shapes:
            shape = (shapes[0] if len(shapes) == 1 else Compound.makeCompound(shapes)).wrapped
        rgb = tuple(self.color.toTuple()[:3]) if self.color is not None else None
        return (self.name, shape, [list(r) for r in self.loc.rows], rgb, [ch._step_tree() for ch in self.children])

    def toCompound(self) -> Compound:
        shapes = list(self.shapes)
        shapes.extend(child.toCompound() for child in self.children)
        return Compound.makeCompound(shapes).locate(self.loc)

    def export(self, path: str, exportType: Optional[str] = None, mode: str = "default", tolerance: float = 0.1, angularTolerance: float = 0.1, **kwargs) -> "Assembly":
        from . import exporters

        exporters.export_assembly(self, path, exportType, tolerance=tolerance, angularTolerance=angularTolerance)
        return self

    save = export

    @classmethod
    def importStep(cls, path: str) -> "Assembly":
        from . import importers

        return importers.importStepAssembly(path)

    @classmethod
    def load(cls, path: str, importType: Optional[str] = None) -> "Assembly":
        return cls.importStep(path)
