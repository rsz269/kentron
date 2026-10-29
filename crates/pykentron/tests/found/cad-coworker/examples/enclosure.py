"""シンプルなネジ穴付き蓋なし筐体（cad-coworker の少し複雑な動作確認用）。

- 7 パラメータをグループ分け（外形 / 取り付け）
- 角の R + 内側くり抜き + 底面四隅のネジ穴
"""

import cadquery as cq

W = 80.0          # 幅 [mm] @min=40 @max=150 @step=1 @group=外形
D = 60.0          # 奥行き [mm] @min=40 @max=150 @step=1 @group=外形
H = 30.0          # 高さ [mm] @min=15 @max=80 @step=1 @group=外形
WALL = 2.0        # 壁厚 [mm] @min=1 @max=5 @step=0.5 @group=外形
FILLET = 3.0      # 角R [mm] @min=0 @max=10 @step=0.5 @group=外形
HOLE = 3.5        # ネジ穴径 [mm] @min=2 @max=6 @step=0.5 @group=取り付け
HOLE_INSET = 5.0  # 穴位置（縁から） [mm] @min=2 @max=20 @step=0.5 @group=取り付け

# 外形（角R付き）
outer = cq.Workplane("XY").box(W, D, H).edges("|Z").fillet(FILLET)

# 内側をくり抜く（壁厚 WALL の凹みを作る）
inner = (
    cq.Workplane("XY")
    .workplane(offset=WALL)
    .box(W - 2 * WALL, D - 2 * WALL, H, centered=(True, True, False))
)
result = outer.cut(inner)

# 底面の四隅にネジ穴を貫通
hole_x = W / 2 - HOLE_INSET
hole_y = D / 2 - HOLE_INSET
result = (
    result.faces("<Z")
    .workplane()
    .pushPoints(
        [(hole_x, hole_y), (-hole_x, hole_y), (hole_x, -hole_y), (-hole_x, -hole_y)]
    )
    .hole(HOLE)
)

show_object(result)  # noqa: F821 (CQGI が注入する)
