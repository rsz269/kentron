"""Cross-section through the lip and diffuser hook at the side of the collar -> out/<size>/snap_section.png
Run: uv run --no-project --with cadquery --with matplotlib snap_section.py"""
import cadquery as cq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import build

slab = cq.Workplane("XY").box(1000, 0.2, 1000)
x0 = build.AO - 22
parts = {
    "skin": (build.skin(), "#ff7a59"),
    "lip (TPU)": (build.lip(), "#f2c230"),
    "diffuser cover": (build.diffuser(), "#d8d2a0"),
}
fig, ax = plt.subplots(figsize=(9, 6), facecolor="#14141a")
ax.set_facecolor("#14141a")
for name, (shape, color) in parts.items():
    cut = shape.intersect(slab)
    if not cut.vals():
        continue
    v, t = cut.val().tessellate(0.05)
    p = np.array([[q.x, q.z] for q in v])
    for tri in t:
        ax.fill(p[list(tri), 0], p[list(tri), 1], color=color, linewidth=0.2, edgecolor=color)
    ax.plot([], [], color=color, label=name, linewidth=6)
ax.set_xlim(x0, build.AO + 6); ax.set_ylim(-48, -16)
ax.set_aspect("equal"); ax.axis("off")
ax.legend(facecolor="#1f1f29", labelcolor="#e8e8ee", loc="upper left", fontsize=9)
ax.set_title("section at the side of the collar: bead on the cover's hook locks into the groove in the lip",
             color="#e8e8ee", fontsize=10, family="monospace")
fig.savefig(build.OUT / "snap_section.png", dpi=130, facecolor=fig.get_facecolor(), bbox_inches="tight")
print("wrote", build.OUT / "snap_section.png")
