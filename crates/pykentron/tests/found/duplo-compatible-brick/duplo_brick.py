"""
DUPLO互換ブロック生成スクリプト (CadQuery)
寸法は実物の実測値ベース:
  外形 31.69mm / ポッチ外径 9.38mm / ポッチ肉厚 1.79mm
  外壁 1.52mm / 裏の筒 内径 10.69mm / 外壁外面〜リブ奥 8.43mm

使い方:
  uv run --python 3.12 --with cadquery python duplo_brick.py   # out/ に STL・STEP を書き出す
  別のサイズは brick(nx, ny) を呼ぶ（例: brick(2, 3)）。
嵌合（STUD_OFFSET / CLEARANCE）と縮み（XY_SHRINK）は Bambu Lab A2L + PETG で合わせた値。
ほかの機種・素材では、高さ半分の試し片を値違いで刷って本物に挿し、合わせ直す。

LEGO / DUPLO は LEGO Group の商標。このスクリプトは非公式で、LEGO Group とは関係ない。
"""
import cadquery as cq
import math

# ---- 基本寸法 (実物) ----
PITCH       = 16.0     # ポッチ間隔
GAP         = 0.155    # 片側の逃げ (16*2 - 31.69)/2
HEIGHT      = 19.2     # ブロック高さ (ポッチ除く)
WALL        = 1.6      # 外壁 (実測1.52 → FDMで0.4mmノズル4周に合わせ1.6)
TOP         = 2.0      # 天板厚
STUD_D      = 9.38     # ポッチ外径
STUD_WALL   = 1.79     # ポッチ肉厚
STUD_H      = 4.5      # ポッチ高さ
TUBE_ID     = 10.69    # 裏の筒 内径
RIB_W       = 1.2      # 内壁リブ幅
EDGE_R      = 1.0      # 外周の縦エッジ丸め
TOP_R       = 0.6      # 天板エッジ丸め
LABEL_SIZE  = 5.0      # 側面の刻印の文字高さ
LABEL_DEPTH = 0.4      # 刻印の深さ（外壁 1.6 のうち）

# ---- 嵌合の既定 (A2L / PETG で 2026-10-04 の t 系テストから決めた) ----
CLEARANCE   = -0.15    # 裏（筒・リブ）を 0.15 きつく。t0.15 が下側ちょうどよかった
STUD_OFFSET = 0.30     # ポッチを 0.30 太く。t0.30 が上側ちょうどよかった
# XY の縮み。設計 31.69 の 2x2 が PETG で 31.51〜31.52 に出た（2026-10-05 実測）→ 1 − 31.515/31.69。
# 外形とポッチ・筒の位置（ピッチ）だけを 1/(1−XY_SHRINK) 倍して、縮んだ後に実物の寸法に戻す。
# ポッチ・筒の径は縮みを含んだ実測で STUD_OFFSET / CLEARANCE を決めたので、倍率をかけない
XY_SHRINK   = 0.0055


def brick(nx=2, ny=2, height=HEIGHT, clearance=CLEARANCE, stud_offset=STUD_OFFSET, label=None):
    """
    clearance  : 裏側の嵌合をゆるくする量(mm)。きつい時は +0.05〜0.2、ゆるい時はマイナス
    stud_offset: ポッチ径の補正(mm)。プリンタで太る時はマイナス
    label      : -Y の側面に彫る文字（嵌合テストの見分け用）
    """
    k = 1 / (1 - XY_SHRINK)
    pitch, gap = PITCH * k, GAP * k
    L = nx * pitch - 2 * gap
    W = ny * pitch - 2 * gap
    stud_d = STUD_D + stud_offset

    # 本体（外形）
    body = (cq.Workplane("XY").box(L, W, height, centered=(True, True, False))
            .edges("|Z").fillet(EDGE_R)
            .faces(">Z").edges().fillet(TOP_R))
    # 中空
    cavity = (cq.Workplane("XY")
              .box(L - 2 * WALL, W - 2 * WALL, height - TOP, centered=(True, True, False)))
    body = body.cut(cavity)

    # ポッチ中心座標
    xs = [(i - (nx - 1) / 2) * pitch for i in range(nx)]
    ys = [(j - (ny - 1) / 2) * pitch for j in range(ny)]
    studs = [(x, y) for x in xs for y in ys]

    # ポッチ（中空・先端面取り）
    for (x, y) in studs:
        s = (cq.Workplane("XY").workplane(offset=height).center(x, y)
             .circle(stud_d / 2).circle(stud_d / 2 - STUD_WALL)
             .extrude(STUD_H)
             .faces(">Z").edges("not %LINE").chamfer(0.4))
        body = body.union(s)

    # 裏の筒（4つのポッチの間に入る）
    inner_h = height - TOP
    tube_od = math.hypot(PITCH, PITCH) - STUD_D - 2 * clearance   # 径なので倍率をかけない
    for i in range(nx - 1):
        for j in range(ny - 1):
            cx = (xs[i] + xs[i + 1]) / 2
            cy = (ys[j] + ys[j + 1]) / 2
            t = (cq.Workplane("XY").center(cx, cy)
                 .circle(tube_od / 2).circle(TUBE_ID / 2)
                 .extrude(inner_h + 0.01))
            body = body.union(t)

    # 内壁リブ（ポッチの外側を押さえる）
    # リブ先端 = ポッチの外接線 + clearance
    stud_edge_from_wall = (pitch / 2 - STUD_D / 2) - gap   # 外面からポッチ外周まで
    rib_depth = stud_edge_from_wall - WALL - clearance
    if rib_depth > 0.2:
        def rib(cx, cy, along_x, sign):
            if along_x:   # 長辺(±Y壁)のリブ
                w = (cq.Workplane("XY")
                     .center(cx, sign * (W / 2 - WALL - rib_depth / 2 + 0.01))
                     .rect(RIB_W, rib_depth + 0.02).extrude(inner_h + 0.01))
            else:
                w = (cq.Workplane("XY")
                     .center(sign * (L / 2 - WALL - rib_depth / 2 + 0.01), cy)
                     .rect(rib_depth + 0.02, RIB_W).extrude(inner_h + 0.01))
            return w
        for x in xs:
            body = body.union(rib(x, 0, True, 1)).union(rib(x, 0, True, -1))
        for y in ys:
            body = body.union(rib(0, y, False, 1)).union(rib(0, y, False, -1))

    # 側面の刻印。XZ 面の法線は -Y なので、負の距離で壁の中へ彫る
    if label:
        mark = (cq.Workplane("XZ", origin=(0, -W / 2, 0)).center(0, height / 2)
                .text(label, LABEL_SIZE, -LABEL_DEPTH, combine=False))
        body = body.cut(mark)

    return body


if __name__ == "__main__":
    import os
    out = os.path.join(os.path.dirname(__file__), "out")
    os.makedirs(out, exist_ok=True)

    b = brick(2, 2)
    cq.exporters.export(b, f"{out}/duplo_2x2.stl", tolerance=0.01, angularTolerance=0.1)
    cq.exporters.export(b, f"{out}/duplo_2x2.step")

    b = brick(2, 4)
    cq.exporters.export(b, f"{out}/duplo_2x4.stl", tolerance=0.01, angularTolerance=0.1)

    # 確認片（ハーフ高さで時短）: 既定の嵌合 + XY の縮み補正。2x4 はポッチ 4 列のピッチが合うかを見る。
    # 刻印 "s55" = XY_SHRINK 0.55% の補正入り（補正前の "30/15" と見分ける）
    for nx, ny in ((2, 2), (2, 4)):
        b = brick(nx, ny, height=9.6, label="s55")
        cq.exporters.export(b, f"{out}/fit_test_{nx}x{ny}_half_s0.55.stl",
                            tolerance=0.01, angularTolerance=0.1)
    print("done")
