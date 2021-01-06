from CGAL.CGAL_Polyhedron_3 import Polyhedron_3
from CGAL.CGAL_Kernel import Point_3
from CGAL.CGAL_Polygon_mesh_processing import Side_of_triangle_mesh
from CGAL.CGAL_Kernel import ON_BOUNDED_SIDE, ON_BOUNDARY, ON_UNBOUNDED_SIDE
import math
import os

"""
case: 'cartesian', 'hexdens', 'cylindric'
cylindric: cylindric distribution along x-axis
"""


class PlacePoints:
    def __init__(self, filename, case, radius):
        self.radius = radius
        self.case = case
        # build polyhedron from file
        self.file = os.getcwd() + filename
        self.polyhedron = Polyhedron_3(self.file)
        assert self.polyhedron.is_closed()

        self.sideof = Side_of_triangle_mesh(self.polyhedron)

        # coordinate ranges x,y,z or x,r,phi
        self.range_i = None
        self.range_j = None
        self.range_k = None
        self.range()

        # Points inside
        self.points = []

        # define points
        if case == 'cartesian':
            self.cartesian()
        elif case == 'hexdens':
            self.hexagonal_densest()
        elif case == 'cylindric':
            self.cylindric()
        else:
            raise ValueError('Unknown case "', case, '"')

    def range(self):
        if self.case == 'cartesian' or self.case == 'hexdens':
            x_range, y_range, z_range = [0., 0.], [0., 0.], [0., 0.]
            for p in self.polyhedron.points():
                x_range[0], x_range[1] = min(x_range[0], p.x()), max(x_range[1], p.x())
                y_range[0], y_range[1] = min(y_range[0], p.y()), max(y_range[1], p.y())
                z_range[0], z_range[1] = min(z_range[0], p.z()), max(z_range[1], p.z())
            self.range_i = x_range
            self.range_j = y_range
            self.range_k = z_range
        elif self.case == 'cylindric':
            x_range, r_range, phi_range = [0., 0.], [0., 0.], [0., 2*math.pi]
            for p in self.polyhedron.points():
                x_range[0], x_range[1] = min(x_range[0], p.x()), max(x_range[1], p.x())
                r_range[1] = max(r_range[1], math.sqrt(p.y()**2+p.z()**2))
            self.range_i = x_range
            self.range_j = r_range
            self.range_k = phi_range
            # axis + transformation matrix

    def cartesian(self):
        print('\nDistribution: cartesian')
        print('Radius particle: ', self.radius)
        dx = 2*self.radius  # dx equals distance particles (2*radius)
        nb_points_i = math.ceil(range_diff(self.range_i) / dx) + 1
        nb_points_j = math.ceil(range_diff(self.range_j) / dx) + 1
        nb_points_k = math.ceil(range_diff(self.range_k) / dx) + 1

        print('Number of points to check: ', nb_points_i*nb_points_j*nb_points_k)

        for i in range(nb_points_i):
            for j in range(nb_points_j):
                for k in range(nb_points_k):
                    px = self.range_i[0] + i * dx
                    py = self.range_j[0] + j * dx
                    pz = self.range_k[0] + k * dx
                    p = Point_3(px, py, pz)
                    if self.sideof.bounded_side(p) == ON_BOUNDED_SIDE or self.sideof.bounded_side(p) == ON_BOUNDARY:
                        self.points.append(p)
        print('Number of points inside polyhedron: ', len(self.points))

    def cylindric(self):
        print('\nDistribution: cylindric')
        print('Radius particle: ', self.radius)
        dx = 2 * self.radius
        points_to_check = 0
        nb_points_i = math.ceil(range_diff(self.range_i) / dx) + 1  # dx
        nb_points_j = math.ceil(range_diff(self.range_j) / dx) + 1  # dr

        # calc points_to_check
        for i in range(nb_points_i):
            for j in range(nb_points_j):
                if j != 0:
                    radius = j*dx
                    dphi = dx/radius
                    nb_points_k = math.floor(2*math.pi/dphi) + 1
                else:
                    nb_points_k = 1
                points_to_check += nb_points_k
        print('\nNumber of points to check: ', points_to_check)

        # check points
        for i in range(nb_points_i):
            for j in range(nb_points_j):
                if j != 0:
                    radius = j*dx
                    dphi = dx/radius
                    nb_points_k = math.floor(2*math.pi/dphi) + 1
                    dphi = 2*math.pi/nb_points_k
                else:
                    dphi = 0
                    nb_points_k = 1
                for k in range(nb_points_k):
                    px = self.range_i[0] + i * dx
                    py = math.cos(k*dphi) * j * dx
                    pz = math.sin(k*dphi) * j * dx
                    p = Point_3(px, py, pz)
                    if self.sideof.bounded_side(p) == ON_BOUNDED_SIDE or self.sideof.bounded_side(p) == ON_BOUNDARY:
                        self.points.append(p)
        print('Number of points inside polyhedron: ', len(self.points))

    def hexagonal_densest(self):
        print('\nDistribution: hexagonal densest package')
        print('Radius particle: ', self.radius)
        dx = self.radius  # dx equals radius of particle
        nb_points_i = math.ceil(range_diff(self.range_i) / (2*dx)) + 1
        nb_points_j = math.ceil(range_diff(self.range_j) / (math.sqrt(3)*dx)) + 1
        nb_points_k = math.ceil(range_diff(self.range_k) / (2*math.sqrt(6)*dx/3)) + 1
        print('\nNumber of points to check: ', nb_points_i * nb_points_j * nb_points_k)

        for i in range(nb_points_i):
            for j in range(nb_points_j):
                for k in range(nb_points_k):
                    px = self.range_i[0] + (2*i + ((j+k) % 2)) * dx
                    py = self.range_j[0] + (math.sqrt(3)*(j+1/3*(k % 2))) * dx
                    pz = self.range_k[0] + (2*math.sqrt(6)*k/3) * dx
                    p = Point_3(px, py, pz)
                    if self.sideof.bounded_side(p) == ON_BOUNDED_SIDE or self.sideof.bounded_side(p) == ON_BOUNDARY:
                        self.points.append(p)
        print('Number of points inside polyhedron: ', len(self.points))


def range_diff(r):
    return abs(r[1]-r[0])
