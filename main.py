from collections import namedtuple

from curves import Pallas
from fields import Fq, Fp
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
    crs1 = setup(8)
    poly1 = Poly([Fq(x) for x in [1, 2, 3, 4, 5, 6, 7, 8]])  # `a` is the coeffs of poly
    x = Fq(1234567890)
    v = poly1.eval(x)
    r1 = Fq.rnd()
    poly_commit = commit(crs1, poly1.coeffs, r1)

    # verifier sends random element U
    u = Pallas.rnd()

    # both parties calculate
    poly_commit_prime = poly_commit + v * u

    # TODO --> top of page 8  ---> need to polish up range limits

    g_prime = crs1.g
    a_prime = poly1.coeffs
    b_prime = [x**i for i in range(len(poly1.coeffs))]

    L = []; R = []; uu = []
    for j in range(len(crs1.g).bit_length() - 1, 0, -1):
        bound = 2**j // 2
        lj = Fq.rnd(); rj = Fq.rnd()
        xx = Poly.dot(a_prime[0:bound], g_prime[bound:]) + lj * crs1.h + Poly.dot(a_prime[0:bound], b_prime[bound::]) * u
        L.append(xx)
        yy = Poly.dot(a_prime[bound:], g_prime[0:bound]) + rj * crs1.h + Poly.dot(a_prime[bound:], b_prime[0:bound]) * u
        R.append(yy)
        uj = Fq.rnd()  # TODO should be drawn from I (challenge space)
        uu.append(uj)
        a_prime = [ahi / uj + alo * uj for (ahi, alo) in zip(a_prime[bound:], a_prime[0:bound])]
        b_prime = [blo / uj + bhi * uj for (blo, bhi) in zip(b_prime[0:bound], b_prime[bound:])]
        g_prime = [glo * (Fq(1)/uj) + ghi * uj for (glo, ghi) in zip(g_prime[0:bound], g_prime[bound:])]


        print('hello ', j, bound, len(a_prime))

    uu.reverse()
    # works
    s = [Fq(1)] * 8
    #u = [3, 6, 9]
    for i in range(2**len(uu)):
        for j in range(len(uu)):
            #print(i & 2**j)
            if i & 2**j != 0:
                s[i] = s[i] * uu[j]
            else:
                s[i] = s[i] / uu[j]
    print("s: ", s)

    g0 = Poly.dot(s, crs1.g)
    print("g0 ", g0)
    print("gp ", g_prime)

    b0 = Poly.dot(s, [x**i for i in range(len(poly1.coeffs))])
    print("b0 ", b0)
    print("bp ", b_prime)
    # omg!!