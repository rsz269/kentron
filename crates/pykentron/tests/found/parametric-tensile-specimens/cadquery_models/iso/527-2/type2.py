import cadquery as cq
import math
import os

# 1. User Parameters (Type 2) 
params = {
    "type_name": "ISO_527-2_Type_2",
    "overall_length": 200.0,       # L0
    "narrow_length": 80.0,         # L2 (Parallel narrow section)
    "narrow_width": 10.0,          # b2
    "grip_width": 20.0,            # b1
    "thickness": 4.0,              # h (User note: can be 4.0 or 10.0)
    "transition_radius": 60.0,     # r
    
    # Informational
    "gauge_length": 50.0,          # L1
    "grip_length_target": 60.0     # Reference only
}

# 2. Geometric Calculations 
def generate_type_2(p):
    # Extract variables
    L0 = p["overall_length"]
    L2 = p["narrow_length"]
    b_narrow = p["narrow_width"]
    b_grip   = p["grip_width"]
    R  = p["transition_radius"]
    h  = p["thickness"]
    
    # Y-Coordinates
    y_narrow = b_narrow / 2.0  # 5.0 mm
    y_grip   = b_grip / 2.0    # 10.0 mm
    
    # Calculate Transition Geometry (dx)
    dy = y_grip - y_narrow      # 5.0 mm
    # Pythagoras: R^2 = dx^2 + (R - dy)^2
    dx = math.sqrt(R**2 - (R - dy)**2) # ~23.98 mm

    # X-Coordinates
    x_start_arc = L2 / 2.0      # 40.0 mm (End of parallel part)
    x_end_arc   = x_start_arc + dx
    x_total     = L0 / 2.0      # 100.0 mm
    
    # Check Resulting Grip Length
    actual_grip_len = x_total - x_end_arc
    
    print(f"--- TYPE 2 GEOMETRY REPORT ---")
    print(f"Parallel Section (L2): {L2} mm")
    print(f"Transition Length (dx): {dx:.2f} mm")
    print(f"Shoulder Start (x): {x_end_arc:.2f} mm")
    print(f"Remaining Grip Length: {actual_grip_len:.2f} mm")
    
    # 3. Modeling 
    # Sketch Top-Right Quadrant
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

    # Extrude and Union
    quarter_solid = q_sketch.extrude(h)
    half_solid    = quarter_solid.union(quarter_solid.mirror("YZ"))
    final_body    = half_solid.union(half_solid.mirror("XZ"))
    
    return final_body

# 4. Execution 
final_model = generate_type_2(params)

# 5. Render in CQ-Editor 
if 'show_object' in globals():
    show_object(final_model, name=params["type_name"])

# 6. Export (Optional) 
output_folder = r"C:\Users\taner\Downloads"
file_path = os.path.join(output_folder, params["type_name"] + ".step")

cq.exporters.export(final_model, file_path)
