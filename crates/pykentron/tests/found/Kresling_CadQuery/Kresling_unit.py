import cadquery as cq
from math import cos, sin, radians
import numpy as np
from cadquery import exporters

# Design parameters
num_sides = 6 # Number of sides 
diameter = 50  # Outer diameter
height = 35  # Dimensions are in millimeter 
theta = 30  # Twisting angle of the upper base in degrees 
base_extrusion = 1  # Base thickness
hole_diameter = 20  # Diameter of the through-hole
scale_factor_top = 0.6  # Scale factor of the top small triangles relative to the surface triangles 
scale_factor_base = 0.8 # Scale factor for the base of small triangles relative to the surface triangles
triangle_thickness = 0.6
small_triangle_thickness = 0.4


rotation_angle = 360 / num_sides
theta_rad = radians(theta)
radius = diameter / 2
offset_small = triangle_thickness + small_triangle_thickness
rotation_angle_rad = radians(rotation_angle)

# Creation of the bases
base1 = cq.Workplane("XY").polygon(num_sides, diameter).extrude(base_extrusion)
base2 = base1.translate((0, 0, height)).rotate((0, 0, 0), (0, 0, 1), theta)

# Definition of the points
P1 = np.array([radius, 0, base_extrusion])
P2 = np.array([radius * cos(theta_rad), radius * sin(theta_rad), height])
P3 = np.array([radius * cos(rotation_angle_rad - theta_rad), -radius * sin(rotation_angle_rad - theta_rad), height])
P4 = np.array([radius * cos(rotation_angle_rad), -radius * sin(rotation_angle_rad), base_extrusion])

P1_tuple = tuple(P1)
P2_tuple = tuple(P2)
P3_tuple = tuple(P3)
P4_tuple = tuple(P4)

# Function to offset the points outward
def offset_points(points, triangle_thickness):
    offset_points = []
    for i in range(len(points)):
        prev_point = points[i - 1]
        next_point = points[(i + 1) % len(points)]

        # Calculate the edge vectors
        edge1 = next_point - points[i]
        edge2 = prev_point - points[i]

        # Calculate the normal vector by averaging the normals to the edges
        normal = np.cross(edge1, edge2)
        normal = normal / np.linalg.norm(normal)  # Normalize the vector
        offset_points.append(points[i] + normal * triangle_thickness)

    return offset_points

# Function to scale down a triangle
def scale_triangle(triangle, scale_factor_base):
    center = np.mean(triangle, axis=0)
    scaled_triangle = [scale_factor_base * (point - center) + center for point in triangle]
    return scaled_triangle

small_triangle_points = scale_triangle([P1, P2, P4], scale_factor_top)
small_triangle2_points = scale_triangle([P2, P4, P3], scale_factor_top)
scaled_triangle_points = scale_triangle([P1, P2, P4], scale_factor_base)
scaled_triangle22_points = scale_triangle([P2, P4, P3], scale_factor_base)

# Offset the points inward (negative offset distance)
small_triangle_points = offset_points(small_triangle_points, -offset_small)
small_triangle_tuples = [tuple(point) for point in small_triangle_points]
small_triangle2_points = offset_points(small_triangle2_points, -offset_small)
small_triangle2_tuples = [tuple(point) for point in small_triangle2_points]
scaled_triangle2_points = offset_points([P2, P4, P3], -triangle_thickness)
scaled_triangle2_tuples = [tuple(point) for point in scaled_triangle2_points]
scaled_triangle_tuples = [tuple(point) for point in scaled_triangle_points]
scaled_triangle22_tuples = [tuple(point) for point in scaled_triangle22_points]
scaled_triangle1_points = offset_points([P1, P2, P4], -triangle_thickness)
scaled_triangle1_tuples = [tuple(point) for point in scaled_triangle1_points]

# Function to create connecting solids
def create_connecting_solids(triangle1_points, scaled_triangle1_points):
    solid = (
        cq.Workplane("XY")
        .polyline(triangle1_points)
        .close()
        .polyline(scaled_triangle1_points)
        .close()
        .loft(combine=True)
    )
    return solid

# Create connecting solids
small_triangle_in = create_connecting_solids(scaled_triangle_tuples, small_triangle_tuples)
small_triangle_out = create_connecting_solids(scaled_triangle22_tuples, small_triangle2_tuples)

triangle2_solid = create_connecting_solids([P2_tuple, P4_tuple, P3_tuple], scaled_triangle2_tuples)
triangle1_solid = create_connecting_solids([P1_tuple, P2_tuple, P4_tuple], scaled_triangle1_tuples)
hinge = create_connecting_solids([P2_tuple, scaled_triangle1_tuples[-2], scaled_triangle2_tuples[0]], [P4_tuple, scaled_triangle1_tuples[-1], scaled_triangle2_tuples[-2]])

# Applying rotations about the z axis
def apply_rotations(shape, num_rotations, rotation_angle, axis_point, axis_direction):
    shapes = [shape]
    for _ in range(num_rotations - 1):
        shapes.append(shapes[-1].rotate(axis_point, axis_direction, rotation_angle))
    return shapes

small_triangle_in_rotated = apply_rotations(small_triangle_in, num_sides, rotation_angle, (0, 0, 0), (0, 0, 1))
small_triangle_out_rotated = apply_rotations(small_triangle_out, num_sides, rotation_angle, (0, 0, 0), (0, 0, 1))
triangle2_solid_rotated = apply_rotations(triangle2_solid, num_sides, rotation_angle, (0, 0, 0), (0, 0, 1))
triangle1_solid_rotated = apply_rotations(triangle1_solid, num_sides, rotation_angle, (0, 0, 0), (0, 0, 1))
hinge_rotated = apply_rotations(hinge, num_sides, rotation_angle, (0, 0, 0), (0, 0, 1))

# Combining all the objects to form one
combined_object = base1.union(base2)

for solid in small_triangle_in_rotated:
    combined_object = combined_object.union(solid)
for solid2 in small_triangle_out_rotated:
    combined_object = combined_object.union(solid2)
for i in triangle2_solid_rotated:
    combined_object = combined_object.union(i)
for j in triangle1_solid_rotated:
    combined_object = combined_object.union(j)
for k in hinge_rotated:
    combined_object = combined_object.union(k)

# Calculate total height
total_height = height  # if you want the hole to be a through all : total_height = height + 2*base_extrusion

# Create the through-hole
combined_object = combined_object.faces(">Z").workplane().hole(hole_diameter, depth=total_height)

# Exporting the object as STL file
exporters.export(combined_object, 'Kresling_unit.stl') # If you want step file just modify the extension

# Display the results
show_object(combined_object, name='Combined Extruded Triangles with Through-Hole in 3D')


