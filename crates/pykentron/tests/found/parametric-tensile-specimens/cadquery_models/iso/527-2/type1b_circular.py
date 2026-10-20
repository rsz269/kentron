import cadquery as cq
import math
import os

# ==============================================================================
# 1. PARAMETERS
# ==============================================================================
params = {
    # Specimen Geometry (ISO 527-2 Type 1B)
    "overall_length": 150.0,
    "gauge_length": 110.0,       
    "parallel_length": 60.0,     # Length of the narrow parallel section (L1)
    "gauge_width": 10.0,         # Width of the narrow section (b1)
    "tab_width": 20.0,           # Width of the grip section (b2)
    "thickness": 4.0,            # Specimen thickness (h)
    "transition_radius": 60.0,   # Radius of the shoulder (r)
    
    # Perforation Pattern Constraints
    "hole_radius": 1.5,          # Radius of the perforation holes
    "hole_spacing": 3.8,         # Center to center pitch distance
    "perimeter_clearance": 0.8,  # Minimum solid wall thickness to maintain
}

# ==============================================================================
# 2. GEOMETRY GENERATION
# ==============================================================================
def generate_boundary_compliant_specimen(p):
    # --- A. Extract Geometric Definitions ---
    L_total = p["overall_length"]
    W_narrow = p["gauge_width"]
    W_grip = p["tab_width"]
    L_parallel = p["parallel_length"]
    R = p["transition_radius"]
    H = p["thickness"]
    
    # Calculate semi-widths (distance from neutral axis Y=0)
    semi_w_narrow = W_narrow / 2.0
    semi_w_grip = W_grip / 2.0
    
    # Calculate Transition Tangency
    # Determine the horizontal length (dx) required for the arc R 
    # to transition smoothly between the two widths.
    dy = semi_w_grip - semi_w_narrow
    dx = math.sqrt(R**2 - (R - dy)**2)
    
    # Define Critical Longitudinal Coordinates (X-Axis)
    x_parallel_end = L_parallel / 2.0       # End of parallel gauge
    x_grip_start = x_parallel_end + dx      # Start of grip tab
    x_terminus = L_total / 2.0              # Physical end of specimen

    # B. Construct Base Solid Body 
    print("[INFO] Generative Design: Creating base solid geometry...")
    
    # Generate 2D Profile (Top-Right Quadrant -> Mirrored)
    profile_sketch = (
        cq.Workplane("XY")
        .moveTo(x_parallel_end, semi_w_narrow)
        .radiusArc((x_grip_start, semi_w_grip), -R)     # Concave transition fillet
        .lineTo(x_terminus, semi_w_grip)
        .lineTo(x_terminus, -semi_w_grip)
        .lineTo(x_grip_start, -semi_w_grip)
        .radiusArc((x_parallel_end, -semi_w_narrow), -R)
        .lineTo(-x_parallel_end, -semi_w_narrow)
        .radiusArc((-x_grip_start, -semi_w_grip), -R)
        .lineTo(-x_terminus, -semi_w_grip)
        .lineTo(-x_terminus, semi_w_grip)
        .lineTo(-x_grip_start, semi_w_grip)
        .radiusArc((-x_parallel_end, semi_w_narrow), -R)
        .close()
    )
    
    solid_body = profile_sketch.extrude(H)

    # C. Compute Compliant Coordinates (Boundary Check) 
    print("[INFO] Pattern Logic: Computing boundary-compliant coordinates...")
    
    r_hole = p["hole_radius"]
    pitch = p["hole_spacing"]
    clearance = p["perimeter_clearance"]
    
    # Minimum Distance Constraint:
    # A hole center must be at least (Clearance + Radius) away from any edge.
    min_edge_dist = clearance + r_hole
    
    # Define Scanning Grid
    # Create a grid slightly larger than the bounding box
    grid_count_x = int(L_total / pitch) + 4
    grid_count_y = int(W_grip / pitch) + 4
    
    origin_x = - (grid_count_x * pitch) / 2
    origin_y = - (grid_count_y * pitch) / 2
    
    compliant_points = []
    
    for i in range(grid_count_x):
        for j in range(grid_count_y):
            # Current Candidate Point
            px = origin_x + i * pitch
            py = origin_y + j * pitch
            
            abs_px = abs(px)
            abs_py = abs(py)
            
            # CONSTRAINT 1: Longitudinal Limits (End Tabs) 
            # Check if point violates the ends of the specimen
            if abs_px > (x_terminus - min_edge_dist):
                continue 

            # CONSTRAINT 2: Transverse Limits (Side Walls) 
            # Determine the maximum allowable Y for the current X position
            max_transverse_pos = 0.0
            
            # Zone 1: Narrow Gauge Section
            if abs_px <= x_parallel_end:
                max_transverse_pos = semi_w_narrow - min_edge_dist
                
            # Zone 2: Grip Tab Section
            elif abs_px >= x_grip_start:
                max_transverse_pos = semi_w_grip - min_edge_dist
                
            # Zone 3: Transition Region (Fillet)
            # Analytically calculate the curve height at this X
            else:
                # Arc center coordinates relative to quadrant
                arc_center_y = semi_w_narrow + R
                x_relative = abs_px - x_parallel_end
                
                # Circle Equation: y = y_center - sqrt(R^2 - x^2)
                curve_height = arc_center_y - math.sqrt(R**2 - x_relative**2)
                
                max_transverse_pos = curve_height - min_edge_dist

            # Final Compliance Check
            if abs_py <= max_transverse_pos:
                compliant_points.append((px, py))

    # D. Apply Boolean Operations
    if compliant_points:
        print(f"[INFO] Boolean Operation: Subtracting {len(compliant_points)} features...")
        
        # Create machining tool (all cylinders combined)
        machining_tool = (
            cq.Workplane("XY")
            .pushPoints(compliant_points)
            .circle(r_hole)
            .extrude(H * 3)         # Oversize extrusion for clean cut
            .translate((0, 0, -H))  # Center vertically
        )
        
        final_part = solid_body.cut(machining_tool)
        return final_part
    else:
        print("[WARNING] Pattern Generation Failed: No coordinates fit within the defined clearance.")
        return solid_body

# ==============================================================================
# 3. EXPORT AND VISUALIZATION
# ==============================================================================

# Execute Generation
final_model = generate_boundary_compliant_specimen(params)

# Visualization in CQ-Editor
if 'show_object' in globals():
    show_object(final_model, name="ISO_527_Perforated_Specimen")

# File Output Configuration
output_directory = r"C:\Users\taner\Downloads"
filename_stl = "ISO_527_Type1b_Circular_v1.stl"
filename_step = "ISO_527_Type1b_Circular_v1.step"

if os.path.exists(output_directory):
    
    # --- 1. STL Export (Standard for 3D Printing) ---
    file_path_stl = os.path.join(output_directory, filename_stl)
    try:
        # Export Settings:
        # tolerance=0.01: High mesh density
        # angularTolerance=0.1: Smooth curvature for DLP
        cq.exporters.export(final_model, file_path_stl, tolerance=0.01, angularTolerance=0.1)
        
        print("-" * 60)
        print(f"[SUCCESS] STL Export Complete.")
        print(f"File: {file_path_stl}")
    except Exception as e:
        print(f"[ERROR] STL Export Failed: {e}")

    # --- 2. STEP Export (Standard for CAD/Engineering) ---
    file_path_step = os.path.join(output_directory, filename_step)
    try:
        cq.exporters.export(final_model, file_path_step)
        print(f"[SUCCESS] STEP Export Complete.")
        print(f"File: {file_path_step}")
        print("-" * 60)
    except Exception as e:
        print(f"[ERROR] STEP Export Failed: {e}")

else:
    print(f"[ERROR] Directory Not Found: {output_directory}")