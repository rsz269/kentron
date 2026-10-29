# Derived script: not one of the hunted corpus's scripts (see derived.tsv and ../README.md).
#
# Derived from: https://github.com/medicationforall/cqindustry at c583078,
#   example/depot/container.py and src/cqindustry/depot/Container.py, Apache License 2.0
#   (the license's notice follows below). It imports the project's own cqindustry and
#   cadqueryhelper (run it with them installed, as the corpus's other cqindustry scripts).
#
# Change (2026-10-29): the original's Container rounds its wedges with
# faces()[0].edges()[1,3].fillet(...) and then faces()[6].edges()[1,3].fillet(...), picking
# faces and edges by their place in the kernel's lists (after the first fillet OCCT's seventh
# face is the wedge's other end, kentron's one of the new fillet faces), so the original's
# result depends on face order and is left out of the corpus counts. Here the edges OCCT's
# indices pick are picked by position: the side edges of the wedge's end at -x, then of its
# end at +x (OCCT's result is unchanged); the rest of the model (the top fillets, the shell
# cut from the body) stays a comparison.
#
# Copyright James Adams (medicationforall).
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file
# except in compliance with the License. You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the
# License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
# either express or implied. See the License for the specific language governing permissions
# and limitations under the License.
import cadquery as cq
from cqindustry.depot import Container


def rounded_wedge(length, height, width, top_length, top_width, side_fillet, top_fillet):
    length_inset = (length - top_length) / 2
    width_inset = (width - top_width) / 2
    wedge = (
        cq.Workplane("XY")
        .wedge(length, height, width, length_inset, width_inset, length - length_inset, width - width_inset)
        .rotate((1, 0, 0), (0, 0, 0), -90)
    )
    # (The edges OCCT's faces()[0].edges()[1,3] and faces()[6].edges()[1,3] pick.)
    sides = wedge.faces("<X").edges("not(>Z or <Z)").fillet(side_fillet)
    sides = sides.faces(">X").edges("not(>Z or <Z)").fillet(side_fillet)
    return sides.faces("Z").fillet(top_fillet)


class ExplicitContainer(Container):
    def make_outline(self):
        body = rounded_wedge(self.length, self.height, self.width, self.top_length, self.top_width, self.side_fillet, self.top_fillet)
        self.outline = body.translate((0, 0, self.height / 2))

    def make_body(self):
        body = rounded_wedge(self.length, self.height, self.width, self.top_length, self.top_width, self.side_fillet, self.top_fillet)
        self.body = body.translate((0, 0, self.height / 2))

    def make_shell(self):
        m = self.material_width
        height = self.height - m
        shell = rounded_wedge(self.length - 2 * m, height, self.width - 2 * m, self.top_length - 2 * m, self.top_width - 2 * m, self.side_fillet, self.top_fillet)
        self.shell = shell.translate((0, 0, height / 2))


bp_container = ExplicitContainer()

bp_container.length = 126
bp_container.width = 114
bp_container.height = 65

bp_container.top_length = 111
bp_container.top_width = 93
bp_container.top_fillet = 8
bp_container.side_fillet = 9

bp_container.material_width = 1

bp_container.make()

ex_container = bp_container.build()

#show_object(ex_container)
cq.exporters.export(ex_container, "stl/depot_container.stl")
