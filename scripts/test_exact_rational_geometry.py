"""Exact geometry primitives independently tested with known rational areas."""
from fractions import Fraction as Q
import unittest
from verify_fixed_zero_rational import area,clipped,triangles


class ExactGeometryTests(unittest.TestCase):
    def test_contact_and_exact_dyadic_overlap(self):
        square=[(Q(0),Q(0)),(Q(1),Q(0)),(Q(1),Q(1)),(Q(0),Q(1))]
        touching=[(x+1,y) for x,y in square]
        self.assertEqual(area(clipped(square,touching)),0)
        epsilon=Q(1,2**40)
        shifted=[(x+1-epsilon,y) for x,y in square]
        self.assertEqual(area(clipped(square,shifted)),epsilon)

    def test_concave_notch_and_area_partition(self):
        shape=[tuple(map(Q,p)) for p in [(0,0),(2,0),(2,1),(1,1),(1,2),(0,2)]]
        parts=triangles(shape)
        self.assertEqual(sum(map(area,parts)),3)
        notch=[tuple(map(Q,p)) for p in [(1,1),(2,1),(2,2),(1,2)]]
        self.assertEqual(sum(area(clipped(t,notch)) for t in parts),0)

    def test_original_self_intersection_is_rejected(self):
        with self.assertRaises(AssertionError):
            triangles([tuple(map(Q,p)) for p in [(0,0),(2,2),(0,2),(2,0)]])


if __name__=='__main__':unittest.main()
