from CGAL.CGAL_Polyhedron_3 import Polyhedron_3
from CGAL.CGAL_Kernel import Point_3
from CGAL.CGAL_Polygon_mesh_processing import Side_of_triangle_mesh
from CGAL.CGAL_Kernel import ON_BOUNDED_SIDE, ON_BOUNDARY, ON_UNBOUNDED_SIDE
import pandas as pd
import math
import os
import sys

"""
case: 'cartesian', 'hexdens', 'cylindric'
cylindric: cylindric distribution along x-axis
"""


class PlacePoints:
    def __init__(self, filename, case, radius_particle):
        self.radius = radius_particle
        self.case = case
        # build polyhedron from file
        self.file = os.getcwd() + filename
        self.polyhedron = Polyhedron_3(self.file)
        assert self.polyhedron.is_closed(), 'Polyhedron not closed'
        self.sideof = Side_of_triangle_mesh(self.polyhedron)

        # coordinate ranges x,y,z or x,r,phi from off-geometry
        self.range_i = None
        self.range_j = None
        self.range_k = None

        # Points inside
        self.points = []  # store CGAL Point_3

        # define points
        if case == 'cartesian':
            self.range_cart()
            self.cartesian()
        elif case == 'hexdens':
            self.range_cart()
            self.hexagonal_densest()
        elif case == 'cylindric':
            self.range_cyl()
            self.cylindric()
        else:
            raise ValueError('Unknown case "', case, '"')

        self.point_coordinates = self.create_points()  # Dataframe containing all point-coordinates

    def range_cart(self):  # defines cartesian coordinate range of off-geometry
        x_range, y_range, z_range = [0., 0.], [0., 0.], [0., 0.]
        for p in self.polyhedron.points():
            x_range[0], x_range[1] = min(x_range[0], p.x()), max(x_range[1], p.x())
            y_range[0], y_range[1] = min(y_range[0], p.y()), max(y_range[1], p.y())
            z_range[0], z_range[1] = min(z_range[0], p.z()), max(z_range[1], p.z())
        self.range_i = x_range
        self.range_j = y_range
        self.range_k = z_range

    def range_cyl(self):  # defines cylindric coordinate range of off-geometry
        x_range, r_range, phi_range = [0., 0.], [0., 0.], [0., 2 * math.pi]
        for p in self.polyhedron.points():
            x_range[0], x_range[1] = min(x_range[0], p.x()), max(x_range[1], p.x())
            r_range[1] = max(r_range[1], math.sqrt(p.y() ** 2 + p.z() ** 2))
        self.range_i = x_range
        self.range_j = r_range
        self.range_k = phi_range

    def cartesian(self):
        print('\nDistribution: cartesian')
        print('Radius particle: ', self.radius)
        dx = 2*self.radius  # dx equals distance particles (2*radius)
        nb_points_i = math.ceil(range_diff(self.range_i) / dx) + 1
        nb_points_j = math.ceil(range_diff(self.range_j) / dx) + 1
        nb_points_k = math.ceil(range_diff(self.range_k) / dx) + 1

        points_total = nb_points_i*nb_points_j*nb_points_k
        print('Number of points to check: ', points_total)

        progress = 0
        for i in range(nb_points_i):
            for j in range(nb_points_j):
                for k in range(nb_points_k):
                    px = self.range_i[0] + i * dx
                    py = self.range_j[0] + j * dx
                    pz = self.range_k[0] + k * dx
                    p = Point_3(px, py, pz)
                    if self.sideof.bounded_side(p) == ON_BOUNDED_SIDE or self.sideof.bounded_side(p) == ON_BOUNDARY:
                        self.points.append(p)
                    progress += 1
                    if progress % 1000 == 0 or progress == points_total:
                        process(progress, points_total)
        print('\nNumber of points inside polyhedron: ', len(self.points))

    def cylindric(self):
        print('\nDistribution: cylindric')
        print('Radius particle: ', self.radius)
        dx = 2 * self.radius
        points_total = 0
        nb_points_i = math.ceil(range_diff(self.range_i) / dx) + 1  # dx
        nb_points_j = math.ceil(range_diff(self.range_j) / dx) + 1  # dr
        # nb_points_k depends on radius (index j)

        # count points_to_check
        for i in range(nb_points_i):
            for j in range(nb_points_j):
                if j != 0:
                    radius = j*dx
                    dphi = dx/radius
                    nb_points_k = math.floor(2*math.pi/dphi) + 1
                else:
                    nb_points_k = 1
                points_total += nb_points_k
        print('\nNumber of points to check: ', points_total)

        # check points
        progress = 0
        for i in range(nb_points_i):
            for j in range(nb_points_j):
                if j != 0:
                    radius = j*dx
                    dphi = dx/radius  # exact dphi to place two points with distance dx
                    nb_points_k = math.floor(2*math.pi/dphi) + 1
                    dphi = 2*math.pi/nb_points_k  # rounded dphi to fit integer nb_points on circumference
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
                    progress += 1
                    if progress % 1000 == 0 or progress == points_total:
                        process(progress, points_total)
        print('\nNumber of points inside polyhedron: ', len(self.points))

    def hexagonal_densest(self):
        print('\nDistribution: hexagonal densest package')
        print('Radius particle: ', self.radius)
        dx = self.radius  # dx equals radius of particle
        nb_points_i = math.ceil(range_diff(self.range_i) / (2*dx)) + 1
        nb_points_j = math.ceil(range_diff(self.range_j) / (math.sqrt(3)*dx)) + 1
        nb_points_k = math.ceil(range_diff(self.range_k) / (2*math.sqrt(6)*dx/3)) + 1
        points_total = nb_points_i * nb_points_j * nb_points_k
        print('\nNumber of points to check: ', points_total)

        progress = 0
        for i in range(nb_points_i):
            for j in range(nb_points_j):
                for k in range(nb_points_k):
                    px = self.range_i[0] + (2*i + ((j+k) % 2)) * dx
                    py = self.range_j[0] + (math.sqrt(3)*(j+1/3*(k % 2))) * dx
                    pz = self.range_k[0] + (2*math.sqrt(6)*k/3) * dx
                    p = Point_3(px, py, pz)
                    if self.sideof.bounded_side(p) == ON_BOUNDED_SIDE or self.sideof.bounded_side(p) == ON_BOUNDARY:
                        self.points.append(p)
                    progress += 1
                    if progress % 1000 == 0 or progress == points_total:
                        process(progress, points_total)
        print('\nNumber of points inside polyhedron: ', len(self.points))

    def create_points(self):
        out = [[p.x(), p.y(), p.z()] for p in self.points]
        out = pd.DataFrame(out, columns=['x', 'y', 'z'])
        return out

    def write_csv(self, filename):
        self.point_coordinates.to_csv(filename)


def range_diff(r):
    return abs(r[1]-r[0])


def process(points_comp, points_total):
    out = '\r{}% --- {} points of {} checked'.format(int(points_comp/points_total*100), points_comp, points_total)
    sys.stdout.write(out)
    sys.stdout.flush()
