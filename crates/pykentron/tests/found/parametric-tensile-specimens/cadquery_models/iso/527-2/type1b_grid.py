import cadquery as cq
import math
import os

# 1. User Parameters 
params = {
    "overall_length": 150.0,
    "gauge_length": 110.0,
    "parallel_length": 60.0,
    "gauge_width": 10.0,
    "tab_width": 20.0,
    "thickness": 4.0,
    "transition_radius": 60.0,
    
    # Manufacturing Parameters 
    "grid_cell_size": 4.0,       # Size of one repeating unit
    "grid_wall_thickness": 0.8,  # Thickness of the internal lattice ribs
    "perimeter_wall": 0.8,       # Thickness of the solid outer shell
}

def generate_iso_with_wall(p):
    # 1. Extract Dimensions
    L_total = p["overall_length"]
    W_narrow = p["gauge_width"]
    W_grip = p["tab_width"]
    L_parallel = p["parallel_length"]
    R = p["transition_radius"]
    H = p["thickness"]
    
    # 2. Geometry Math
    y_narrow = W_narrow / 2.0
    y_grip = W_grip / 2.0
    dy = y_grip - y_narrow
    dx = math.sqrt(R**2 - (R - dy)**2)
    x_start_arc = L_parallel / 2.0
    x_end_arc = x_start_arc + dx
    x_total = L_total / 2.0

    # 3. Create the 2D Profile (Full Closed Loop)
    q_edge = (
        cq.Workplane("XY")
        .moveTo(0, 0)
        .lineTo(0, y_narrow)
        .lineTo(x_start_arc, y_narrow)
        .radiusArc((x_end_arc, y_grip), -R)
        .lineTo(x_total, y_grip)
        .lineTo(x_total, 0)
        .close()
    )
    
    # Create the full solid profile (Outer Shell)
    quarter_solid = q_edge.extrude(H)
    half_solid = quarter_solid.union(quarter_solid.mirror("YZ"))
    main_body = half_solid.union(half_solid.mirror("XZ"))
    
    # 4. Create the Inner Core (For the Grid) 
    full_profile = (
        cq.Workplane("XY")
        .moveTo(x_start_arc, y_narrow)
        .radiusArc((x_end_arc, y_grip), -R)
        .lineTo(x_total, y_grip)
        .lineTo(x_total, -y_grip)
        .lineTo(x_end_arc, -y_grip)
        .radiusArc((x_start_arc, -y_narrow), -R)
        .lineTo(-x_start_arc, -y_narrow)
        .radiusArc((-x_end_arc, -y_grip), -R)
        .lineTo(-x_total, -y_grip)
        .lineTo(-x_total, y_grip)
        .lineTo(-x_end_arc, y_grip)
        .radiusArc((-x_start_arc, y_narrow), -R)
        .close()
    )
    
    # Shrink the profile by the wall thickness
    inner_core_solid = (
        full_profile
        .offset2D(-p["perimeter_wall"]) 
        .extrude(H)
    )

    # 5. Generate Grid Cutter 
    cell_size = p["grid_cell_size"]
    hole_size = cell_size - p["grid_wall_thickness"]
    
    count_x = int(L_total / cell_size) + 2
    count_y = int(W_grip / cell_size) + 2
    
    raw_grid = (
        cq.Workplane("XY")
        .rarray(cell_size, cell_size, count_x, count_y, center=True)
        .rect(hole_size, hole_size)
        .extrude(H * 2)             # Extrude tall
        .translate((0, 0, -H/2))    # Center in Z
    )
    
    # 6. Boolean Logic 
    constrained_holes = raw_grid.intersect(inner_core_solid)
    final_part = main_body.cut(constrained_holes)
    
    return final_part

# Execution 
final_model = generate_iso_with_wall(params)

# Render 
if 'show_object' in globals():
    show_object(final_model, name="ISO_527_Walled_Grid")

# EXPORT SECTION,
output_folder = r"C:\Users\taner\Downloads"

if not os.path.exists(output_folder):
    print(f"Warning: Directory '{output_folder}' does not exist.")
else:
    # 1. STL Export
    stl_name = "ISO_527_Walled_Grid.stl"
    stl_path = os.path.join(output_folder, stl_name)
    try:
        cq.exporters.export(final_model, stl_path, tolerance=0.01, angularTolerance=0.1)
        print(f"SUCCESS: STL Saved -> {stl_path}")
    except Exception as e:
        print(f"STL Error: {e}")

# STEP Export (Optional)
# Uncomment the lines below carefully
#    step_name = "ISO_527_Walled_Grid.step"   # <--- Aligned with stl_name above
#    step_path = os.path.join(output_folder, step_name)
#    try:
#        cq.exporters.export(final_model, step_path)
#        print(f"SUCCESS: STEP Saved -> {step_path}")
#    except Exception as e:
#        print(f"STEP Error: {e}")
