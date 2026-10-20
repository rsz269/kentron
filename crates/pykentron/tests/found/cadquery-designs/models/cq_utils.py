# BSD 3-Clause License
#
# Copyright (c) 2025, Viktor Samokhin (wowyupiyo@gmail.com)
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
################################################################################
#
# v0.0.2
# Utility module to share code snippets across the designs

import os

import cadquery as cq
from cadquery import Workplane
from ocp_vscode import show_object


def show_models(**kwargs: Workplane):
    """Show models, passed in as named argumets, with OCP Viewer"""
    for name, model in kwargs.items():
        nice_name = name.replace("_", " ").replace("-", " ").title()
        show_object(model, name=nice_name)

def export_models(stl_export=True, step_export=False, **kwargs: Workplane):
    """Export models, passed in as named arguments, to STL and/or STEP formats"""

    if not (stl_export or step_export): return

    # Create export directory
    export_dir = "exports"
    if not os.path.exists(export_dir):
        os.makedirs(export_dir)

    # Use argument names as file names
    for name, model in kwargs.items():
        if stl_export:
            # Export STL
            stl_path = os.path.join(export_dir, f"{name}.stl")
            cq.exporters.export(model, stl_path)
            print(f"Exported {stl_path}")

        if step_export:
            # Export STEP
            step_path = os.path.join(export_dir, f"{name}.step")
            cq.exporters.export(model, step_path)
            print(f"Exported {step_path}")
