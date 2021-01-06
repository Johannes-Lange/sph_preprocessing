from CGAL.CGAL_Polyhedron_3 import Polyhedron_3
from CGAL.CGAL_Kernel import Point_3
from CGAL.CGAL_Kernel import ON_BOUNDED_SIDE, ON_BOUNDARY
from CGAL.CGAL_Polygon_mesh_processing import Side_of_triangle_mesh
import pandas as pd
import os

# from CGAL.CGAL_Polygon_mesh_processing import orient_polygon_soup

write_file = False
datafile = os.getcwd() + '/bunny.off'

# Create input polyhedron
polyhedron = Polyhedron_3(datafile)
print('\nPolyhedron closed: ', polyhedron.is_closed())

# Find range coordinates (x_max, y_max, z_max)
x_range, y_range, z_range = [0., 0.], [0., 0.], [0., 0.]
for p in polyhedron.points():
    x_range[0], x_range[1] = min(x_range[0], p.x()), max(x_range[1], p.x())
    y_range[0], y_range[1] = min(y_range[0], p.y()), max(y_range[1], p.y())
    z_range[0], z_range[1] = min(z_range[0], p.z()), max(z_range[1], p.z())
print('Coordinates ranges:')
print('x:', x_range, '\ny:', y_range, '\nz:', z_range)

# Define points in range
nb_points_direction = 30
print('Number of points to check: ', nb_points_direction ** 3)

(dx, dy, dz) = (abs(x_range[1] - x_range[0]) / (nb_points_direction - 1),
                abs(y_range[1] - y_range[0]) / (nb_points_direction - 1),
                abs(z_range[1] - z_range[0]) / (nb_points_direction - 1))

"""
dx = dy = dz = 0.03
nb_points_dir = math.ceil(range/dx)+1
"""
# Create points and check if inside polyhedron
inst = Side_of_triangle_mesh(polyhedron)

points = []
for x in range(nb_points_direction):
    for y in range(nb_points_direction):
        for z in range(nb_points_direction):
            px = x_range[0] + x * dx
            py = y_range[0] + y * dy
            pz = z_range[0] + z * dz
            p = Point_3(px, py, pz)
            if inst.bounded_side(p) == ON_BOUNDED_SIDE or inst.bounded_side(p) == ON_BOUNDARY:
                points.append(p)

print('Number of points inside polyhedron: ', len(points))

# Write points to csv
out = [[p.x(), p.y(), p.z()] for p in points]
out = pd.DataFrame(out, columns=['x', 'y', 'z'])
if write_file:
    out.to_csv('points_inside.csv')
    print('\nFile written!')
