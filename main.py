import secrets
from collections import namedtuple

from curves import Pallas
from fields import Fp, Fq
from poly import Poly

Crs = namedtuple("Crs", "g h")


def setup(size):
    g = []
    for i in range(size):
        token_g = secrets.token_bytes(16)
        g.append(Pallas.hash_to_curve(i.to_bytes(1, byteorder='little')))  # token_g))
    token_h = secrets.token_bytes(16)
    h = Pallas.hash_to_curve(b'\x11')  # token_h)
    return Crs(g=g, h=h)


def commit(crs, poly, r):
    dd = Poly.dot(crs.g, poly)
    return dd + crs.h * r


if __name__ == '__main__':
    crs = setup(4)

    pp = Poly([Fq(x) for x in [1, 1]]) #), 2, 3, 4]])
    qq = Poly([Fq(x) for x in [1, 2]]) #, 8, 9, 10]])
    rr = Fq(1)
    ss = Fq(2)
    aa = Fq(1234)
    bb = Fq(5678)
    com = commit(crs, pp.coeffs, rr)
    print(com)

    left = commit(crs, pp.coeffs, rr) * aa + commit(crs, qq.coeffs, ss) * bb
    right = commit(crs, [aa * x + bb * y for (x, y) in zip(pp.coeffs, qq.coeffs)], aa * rr + bb * ss)
    # ziz = [aa * x + bb * y for (x, y) in zip(pp.coeffs, qq.coeffs)]
    # tt = aa * pp.coeffs[0] + bb * qq.coeffs[0]
    lx = left.x / left.z
    rx = right.x / right.z
    hmm = crs.h * rr
    assert left == right
