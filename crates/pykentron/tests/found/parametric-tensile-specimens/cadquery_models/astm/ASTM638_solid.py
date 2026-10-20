import cadquery as cq
import math
import os

# ==============================================================================
# 1. STANDARD SPECIFICATIONS (ASTM D638)
# ==============================================================================


ASTM_SPECS = {
    "TYPE_I": {
        "L_tot": 165.0,  # Overall Length (min)
        "L_par": 57.0,   # Length of narrow parallel section
        "W_nar": 13.0,   # Width of narrow section
        "W_grip": 19.0,  # Width of grip section
        "Thick": 3.2,    # Standard Thickness
        "Rad": 76.0      # Large transition radius
    },
    "TYPE_II": {
        "L_tot": 183.0, 
        "L_par": 57.0, 
        "W_nar": 6.0,    # Much narrower than Type I
        "W_grip": 19.0,
        "Thick": 3.2,   
        "Rad": 25.0      
    },
    "TYPE_V": {
        # Small specimen for limited material
        "L_tot": 63.5,  
        "L_par": 9.53, 
        "W_nar": 3.18, 
        "W_grip": 9.53,
        "Thick": 3.2,   
        "Rad": 12.7
    }
}

# ------------------------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------------------------
# Change this variable to "TYPE_I", "TYPE_II", or "TYPE_V" to switch models.
SELECTED_TYPE = "TYPE_II"

# Directory where files will be saved
OUTPUT_DIR = r"C:\Users\taner\Downloads"

# ==============================================================================
# 2. GEOMETRY GENERATE
# ==============================================================================
def generate_astm_specimen(spec_name):
    """
    Generates a solid 3D model for the selected ASTM specimen type.
    It calculates the transition curve mathematically to ensure smooth tangency.
    """
    print(f"Starting generation for {spec_name}...")
    
    # Get parameters for the selected type
    p = ASTM_SPECS[spec_name]
    
    # Unpack values
    l0, l2 = p["L_tot"], p["L_par"]
    w_narrow, w_grip = p["W_nar"], p["W_grip"]
    rad, thk = p["Rad"], p["Thick"]
    
    # Calculate half-widths (distance from centerline to edge)
    y_narrow = w_narrow / 2.0
    y_grip = w_grip / 2.0
    
    # --- Geometric Calculations 
    # We need to calculate the horizontal distance (dx) that the fillet radius consumes.
    # dy is the vertical step size between the narrow part and the grip.
    dy = y_grip - y_narrow
    
    # Pythagorean theorem to find dx: R^2 = dx^2 + (R - dy)^2
    # This ensures the arc is perfectly tangent to both straight sections.
    if rad >= dy:
        dx = math.sqrt(rad**2 - (rad - dy)**2)
    else:
        print(f"Warning: Radius ({rad}) is too small for the width change ({dy}).")
        dx = dy # Fallback to linear transition if geometry is invalid

    # Define critical X-axis coordinates relative to the center (0,0)
    x_narrow_end = l2 / 2.0          # Where the straight narrow part ends
    x_arc_end = x_narrow_end + dx    # Where the curve ends and grip starts
    x_total_end = l0 / 2.0           # The physical end of the specimen
    
    # Check if the overall length is long enough to contain the geometry
    if x_total_end < x_arc_end:
        print("Warning: Overall length (L0) is too short! Extending grips automatically.")
        x_total_end = x_arc_end + 5.0 # Add 5mm grip buffer

    # --- Sketching the Profile 
    # We draw the top-right quadrant first.
    edge = (
        cq.Workplane("XY")
        .moveTo(0, 0)
        .lineTo(0, y_narrow)                # Start at center, go up
        .lineTo(x_narrow_end, y_narrow)     # Draw the parallel narrow section
        .radiusArc((x_arc_end, y_grip), -rad) # Draw the concave transition arc
        .lineTo(x_total_end, y_grip)        # Draw the grip section
        .lineTo(x_total_end, 0)             # Close back to centerline
        .close()
    )

    # --- Creating the Solid 
    # 1. Extrude the 2D profile to the specified thickness
    quarter_solid = edge.extrude(thk)
    
    # 2. Mirror across YZ plane (creates left side)
    # 3. Mirror across XZ plane (creates bottom side)
    # We use union=True to merge them into a single watertight solid immediately.
    full_body = (
        quarter_solid
        .mirror("YZ", union=True)
        .mirror("XZ", union=True)
    )

    return full_body

# ==============================================================================
# 3. EXECUTION AND EXPORT
# ==============================================================================
try:
    # Check if output directory exists, create if not
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # Generate the model
    model = generate_astm_specimen(SELECTED_TYPE)

    # Show in CadQuery Editor (if running inside it)
    if 'show_object' in globals():
        show_object(model, name=SELECTED_TYPE)

    # Define file paths
    stl_path = os.path.join(OUTPUT_DIR, f"ASTM_D638_{SELECTED_TYPE}.stl")
    step_path = os.path.join(OUTPUT_DIR, f"ASTM_D638_{SELECTED_TYPE}.step")
    
    # Export to STL
    # Ensure the smooth curves
    cq.exporters.export(model, stl_path, tolerance=0.01, angularTolerance=0.05)
    print(f"Success: STL saved to -> {stl_path}")
    
    # Export to STEP 
    cq.exporters.export(model, step_path)
    print(f"Success: STEP saved to -> {step_path}")

except Exception as e:
    # Catch any errors and print them clearly

    print(f"An error occurred: {e}")
