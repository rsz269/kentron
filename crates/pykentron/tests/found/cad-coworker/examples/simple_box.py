"""シンプルな直方体（cad-coworker MVP の動作確認用）。

cad-coworker が読み込む CadQuery スクリプトの最小サンプル。
各パラメータに `@min @max @step @group` を付けてパーサ仕様を網羅する。
"""

import cadquery as cq

WIDTH = 100.0   # 幅 [mm] @min=20 @max=200 @step=0.5 @group=サイズ
DEPTH = 80.0    # 奥行き [mm] @min=20 @max=200 @step=0.5 @group=サイズ
HEIGHT = 30.0   # 高さ [mm] @min=5 @max=100 @step=0.5 @group=サイズ

result = cq.Workplane("XY").box(WIDTH, DEPTH, HEIGHT)

show_object(result)  # noqa: F821 (CQGI が注入する)
