import cadquery as cq
import math
import os

# 1. User Parameters
params = {
    "overall_length": 150.0,      # mm (minimum)
    "gauge_length": 110.0,        # mm (Distance between grips/shoulders L)
    "parallel_length": 60.0,      # mm (Length of narrow parallel section)
    "gauge_width": 10.0,          # mm (Width of narrow section)
    "tab_width": 20.0,            # mm (Width of grip section)
    "thickness": 4.0,             # mm
    "transition_radius": 60.0,    # mm
    "gauge_mark_spacing": 50.0,   # mm (Info only)
}

# 2. Geometric Calculations
def generate_from_params(p):

    L_total = p["overall_length"]
    W_narrow = p["gauge_width"]
    W_grip = p["tab_width"]
    L_parallel = p["parallel_length"]
    R = p["transition_radius"]
    H = p["thickness"]

    # Y-Coordinates (Heights)
    y_narrow = W_narrow / 2.0
    y_grip = W_grip / 2.0

    # Calculate Transition Length (dx)
    # Pythagoras: R^2 = dx^2 + (R - dy)^2
    dy = y_grip - y_narrow
    dx = math.sqrt(R**2 - (R - dy)**2)

    # X-Coordinates
    # 1. End of the straight narrow part
    x_start_arc = L_parallel / 2.0
    
    # 2. End of the curve (Start of the wide grip)
    x_end_arc = x_start_arc + dx
    
    # 3. Total half-length
    x_total = L_total / 2.0

    # Console Info
    shoulder_dist = x_end_arc * 2
    print(f"--- GENERATION INFO ---")
    print(f"Parallel Length: {L_parallel} mm")
    print(f"Shoulder Distance (Calculated): {shoulder_dist:.2f} mm")
    print(f"Grip Distance (Target): {p['gauge_length']} mm")
    print(f"Status: Grip distance is outside the curve? {'YES' if p['gauge_length'] >= shoulder_dist else 'WARNING: Grips overlap curve'}")

    # 3. Modeling ---
    # Top-Right Quadrant Sketch
    q_sketch = (
        cq.Workplane("XY")
        .moveTo(0, 0)
        .lineTo(0, y_narrow)            # Center -> Narrow Width
        .lineTo(x_start_arc, y_narrow)  # Parallel Section
        
        # Negative R (-R) for Inner/Concave Arc
        .radiusArc((x_end_arc, y_grip), -R) 
        
        .lineTo(x_total, y_grip)        # Grip Section
        .lineTo(x_total, 0)             # Close Loop
        .close()
    )

    # Extrude and Union (It prevents split lines)
    quarter_solid = q_sketch.extrude(H)
    half_solid = quarter_solid.union(quarter_solid.mirror("YZ"))
    final_body = half_solid.union(half_solid.mirror("XZ"))
    
    return final_body

# 4. Execution ---
final_model = generate_from_params(params)

# 5. Render in CQ-Editor
if 'show_object' in globals():
    show_object(final_model, name="ISO_527_From_Params")

# 6. Export (Optional) 
# output_folder = r"C:\Users\taner\Downloads" #Change the selected directory
# file_path = os.path.join(output_folder, "ISO_527_Parametric.step")
# cq.exporters.export(final_model, file_path)