from collections import namedtuple

from curves import Pallas
from fields import Fq
from poly import Poly

Crs = namedtuple("Crs", "g h")


def setup(size):
    g = []
    for i in range(size):
        g.append(Pallas.rnd())
    h = Pallas.rnd()
    return Crs(g=g, h=h)


def commit(crs, poly, r):
    dd = Poly.dot(crs.g, poly)
    return dd + crs.h * r


if __name__ == '__main__':
    crs1 = setup(4)

    pp = Poly([Fq(x) for x in [1, 2, 3, 4]])
    qq = Poly([Fq(x) for x in [7, 8, 9, 10]])
    rr = Fq(9876)
    ss = Fq(7654)
    aa = Fq(1234)
    bb = Fq(5678)
    com = commit(crs1, pp.coeffs, rr)
    print(com)

    left = commit(crs1, pp.coeffs, rr) * aa + commit(crs1, qq.coeffs, ss) * bb
    right = commit(crs1, [aa * x + bb * y for (x, y) in zip(pp.coeffs, qq.coeffs)], aa * rr + bb * ss)
    assert left == right

    # -----> TODO: section 3.1

    # inputs to protocol
    crs1 = setup(4)
    poly1 = Poly([Fq(x) for x in [1, 2, 3, 4]])  # `a` is the coeffs of poly
    x = Fq(1234567890)
    v = poly1.eval(x)
    r1 = Fq.rnd()
    poly_commit = commit(crs1, poly1.coeffs, r1)

    # verifier sends random element U
    u = Pallas.rnd()

    # both parties calculate
    poly_commit_prime = poly_commit + v * u





