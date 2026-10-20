# Copyright © 2025 Jeremy Wright
#
# License: BSD
#
# THIS SOFTWARE IS PROVIDED BY JEREMY WRIGHT AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JEREMY WRIGHT BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
# WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import cadquery as cq

# All units are in mm

#
# Parameters
#
ziptie_bar_thickness = 6.0  # How tall the zip-tie rails on the sides are
ziptie_width = 5.0          # Width of the zip-tie band, change if your ziptie is bigger
ziptie_thickness = 2.0      # Thickness of the zip-tie band, change if your ziptie is bigger

# Main shape
plate = cq.Workplane("XY" ).box(155.0, 155.0, 3.0).edges("|Z").fillet(3.0)

# Add the stiffness/ziptie bars
plate = (plate.faces(">Z")
                .workplane()
                .pushPoints([(0.0, -67),(0.0, 67)])
                .rect(124.0, 5.0)
                .extrude(ziptie_bar_thickness))
plate = (plate.faces(">Z").edges("|Y").chamfer(ziptie_bar_thickness - 0.01, ziptie_bar_thickness / 2.0))
plate = (plate.faces(">Z[-2]")
                .workplane()
                .pushPoints([(-67, 0.0),(67, 0.0)])
                .rect(5.0, 124.0)
                .extrude(ziptie_bar_thickness + 0.001))
plate = (plate.faces(">Z").edges("|X").chamfer(ziptie_bar_thickness - 0.01, ziptie_bar_thickness / 2.0))

# Add ziptie cutouts in bars
plate = (plate.faces(">Y[-2]")
                .workplane(centerOption="CenterOfBoundBox")
                .rarray(10.0, 0.0, 12, 1)
                .rect(ziptie_width, ziptie_thickness)
                .cutThruAll())
plate = (plate.faces(">X[-2]")
                .workplane(centerOption="CenterOfBoundBox")
                .rarray(10.0, 0.0, 12, 1)
                .rect(ziptie_width, ziptie_thickness)
                .cutThruAll())

# Plastic reduction areas
plate = (plate.faces("<Z")
                .workplane(centerOption="CenterOfBoundBox")
                .rarray(30.0, 30.0, 4, 4)
                .polygon(6, 24.0, circumscribed=False)
                .cutThruAll())

# Counter-sunk mounting holes
plate = (plate.faces("<Z")
                .workplane(centerOption="CenterOfBoundBox")
                .pushPoints([(-67.5, -67.5), (-67.5, 67.5), (67.5, -67.5), (67.5, 67.5)])
                .circle(5.0 / 2.0)
                .cutBlind(-3.0, taper=-45))

# Uncomment and change the path to export
# cq.exporters.export(plate, '/path/to/model/nimble_ziptie_topplate.stl')

show_object(plate)