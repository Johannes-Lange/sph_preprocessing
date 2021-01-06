from CGAL.CGAL_Kernel import Point_2, Polygon_2, do_intersect  # Polygon_with_holes_2 missing
from CGAL.CGAL_Kernel import ON_BOUNDED_SIDE, ON_BOUNDARY, ON_UNBOUNDED_SIDE


# TODO: wrapper class, method bounded_side
# Check if point inside (with holes)
def inside_shape(outer, holes, point):
    for hole in holes:
        if hole.bounded_side(point) == ON_BOUNDED_SIDE:
            return ON_UNBOUNDED_SIDE
        elif hole.bounded_side(point) == ON_BOUNDARY:
            return ON_BOUNDARY

    if outer.bounded_side(point) == ON_UNBOUNDED_SIDE:
        return ON_UNBOUNDED_SIDE
    if outer.bounded_side(point) == ON_BOUNDED_SIDE:
        return ON_BOUNDED_SIDE
    if outer.bounded_side(point) == ON_BOUNDARY:
        return ON_BOUNDARY


p1 = Point_2(-1, -1)
p2 = Point_2(-1, 1)
p3 = Point_2(1, 1)
p4 = Point_2(1, -1)
outside = Polygon_2([p1, p2, p3, p4])

p5 = Point_2(-.5, -.5)
p6 = Point_2(-.5, .5)
p7 = Point_2(.5, .5)
p8 = Point_2(.5, -.5)
inside = Polygon_2([p5, p6, p7, p8])  # hole

q = Point_2(10, 10)
print(inside_shape(outside, [inside], q) == ON_BOUNDED_SIDE)
