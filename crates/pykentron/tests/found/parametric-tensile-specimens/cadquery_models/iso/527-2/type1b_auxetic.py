import cadquery as cq
import math
import os

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================
# ISO 527-2 Type 1B Dimensions
PARAMS = {
    "L_tot": 150.0,   "L_par": 60.0,    "W_nar": 10.0,
    "W_grip": 20.0,   "Thick": 4.0,     "Rad": 60.0,
}

# Auxetic Pattern Settings
PATTERN = {
    "cell": 4.0,      # Cell size
    "wall": 0.8,      # Wall thickness
    "depth": 0.3,     # Depth of the re-entrant angle
    "margin": 0.5,    # Safety margin from the edge
    "rotate": True    # True = Vertical Bow-tie (90 deg rotation)
}

OUTPUT_DIR = r"C:\Users\taner\Downloads"

# ==============================================================================
# 2. BOUNDARY
# ==============================================================================
def is_inside_boundary(x, y):
    """
    Mathematically checks if a point (x,y) is inside the ISO 527 shape.
    Uses rounding to strictly solve the 'Left Quadrant' missing issue.
    """
    p = PARAMS
    # Effective margin required
    safe_margin = PATTERN["margin"] + (PATTERN["cell"] / 2.0)
    
    # Use absolute values for perfect symmetry (Left = Right)
    ax = round(abs(x), 5)
    ay = round(abs(y), 5)
    
    # 1. Check Total Length limits
    if ax > ((p["L_tot"] / 2.0) - safe_margin):
        return False

    # 2. Check Width limits (Y-Axis)
    max_y = 0.0

    # Region A: Parallel Narrow Section
    if ax <= (p["L_par"] / 2.0):
        max_y = p["W_nar"] / 2.0

    # Region B: Grip & Transition
    else:
        # Calculate transition parameters
        dy = (p["W_grip"] - p["W_nar"]) / 2.0
        dx = math.sqrt(p["Rad"]**2 - (p["Rad"] - dy)**2)
        x_grip_start = (p["L_par"] / 2.0) + dx

        if ax >= x_grip_start:
            # Straight Grip Section
            max_y = p["W_grip"] / 2.0
        else:
            # Transition Curve (Fillet)
            x_local = ax - (p["L_par"] / 2.0)
            if x_local > p["Rad"]: x_local = p["Rad"] # Clamp
            
            # Curve Formula
            drop = p["Rad"] - math.sqrt(p["Rad"]**2 - x_local**2)
            max_y = (p["W_nar"] / 2.0) + drop

    # Final Check
    return ay <= (max_y - safe_margin)

# ==============================================================================
# 3. GENERATION
# ==============================================================================
def generate_ultimate_specimen():
    print("STEP 1: Generating Base Solid...")
    p = PARAMS
    
    # A. Base Geometry 
    yn, yg = p["W_nar"]/2.0, p["W_grip"]/2.0
    dy = yg - yn
    dx = math.sqrt(p["Rad"]**2 - (p["Rad"] - dy)**2)
    xs, xe, xt = p["L_par"]/2.0, p["L_par"]/2.0 + dx, p["L_tot"]/2.0

    # Draw quadrant
    edge = (
        cq.Workplane("XY").moveTo(0,0).lineTo(0, yn)
        .lineTo(xs, yn).radiusArc((xe, yg), -p["Rad"])
        .lineTo(xt, yg).lineTo(xt, 0).close()
    )
    # Extrude and Mirror
    base = edge.extrude(p["Thick"])
    base = base.mirror("YZ", union=True).mirror("XZ", union=True)

    print("STEP 2: Calculating Pattern Coordinates...")
    
    # B. Define Pattern Shape 
    c = PATTERN["cell"]
    w = PATTERN["wall"]
    h = c - w
    d = h / 2.0
    xi = d - (h * PATTERN["depth"])
    
    # Standard Vertices (Horizontal)
    pts_std = [
        (xi, 0), (d, d), (-d, d), 
        (-xi, 0), (-d, -d), (d, -d)
    ]
    
    # Rotate 90 Degrees if requested
    if PATTERN["rotate"]:
        pts_draw = [(y, x) for x, y in pts_std]
    else:
        pts_draw = pts_std

    # C. Brute Force Grid Loop 
    # Scan range: Covers the entire bounding box
    rx = int((p["L_tot"] / 2.0) / c) + 2
    ry = int((p["W_grip"] / 2.0) / c) + 2
    
    # We create a single Workplane
    cutter_sketch = cq.Workplane("XY")
    count = 0

    # Loop from Negative to Positive
    for i in range(-rx, rx + 1):
        for j in range(-ry, ry + 1):
            
            # Center of the cell
            cx = i * c
            cy = j * c
            
            # Mathematical Check
            if is_inside_boundary(cx, cy):
                # Calculate Absolute Coordinates for this cell
                abs_pts = [(vx + cx, vy + cy) for vx, vy in pts_draw]
                
                # Add this polygon to the stack
                cutter_sketch = cutter_sketch.polyline(abs_pts).close()
                count += 1
                
    print(f"STEP 3: Pattern Generated. Total Cells: {count}")

    if count == 0:
        print("ERROR: No cells fit inside. Reduce margin or cell size.")
        return base

    # --- D. Single Boolean Cut ---
    print("STEP 4: Extruding and Cutting...")
    
    # Extrude all accumulated wires at once
    cut_tool = cutter_sketch.extrude(p["Thick"] * 2).translate((0, 0, -p["Thick"]))
    
    # Final Cut
    final_model = base.cut(cut_tool)
    
    return final_model

# ==============================================================================
# 4. EXECUTION & EXPORT
# ==============================================================================
try:
    # 1. Run Generator
    model = generate_ultimate_specimen()
    
    # 2. Show in CQ
    if 'show_object' in globals():
        show_object(model, name="Final_Auxetic_Specimen")

    # 3. Export to STL
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    
    fpath = os.path.join(OUTPUT_DIR, "ISO527_Auxetic_Vertical_Final.stl")
    cq.exporters.export(model, fpath, tolerance=0.01, angularTolerance=0.05)
    print(f"SUCCESS: File saved to {fpath}")

    # ==========================================================================
    # 4. OPTIONAL: STEP EXPORT (Uncomment to enable)
    # ==========================================================================
    # step_path = os.path.join(OUTPUT_DIR, "ISO527_Auxetic_Vertical_Final.step")
    # cq.exporters.export(model, step_path)
    # print(f"SUCCESS: STEP File saved to {step_path}")

except Exception as e:
    import traceback
    traceback.print_exc()
    print(f"CRITICAL ERROR: {e}")