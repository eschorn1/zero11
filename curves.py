from fields import Fp, Fq


class Curve:

    # Constants are to be defined in specific Field subclasses
    b3 = None  # Code will fail if Field is directly utilized

    def __init__(self, x, y, z):
        assert type(x) == Fp or type(x) == Fq
        assert type(x) == type(y) and type(y) == type(z)
        self.x, self.y, self.z = x, y, z

    def __eq__(self, other):
        assert type(self) == type(other)
        if self.z == (type(self.z))(0): return other.z == (type(other.z))(0)
        if other.z == (type(other.z))(0): return False
        return self.x/self.z == other.x/other.z and self.y/self.z == other.y/other.z

    # See https: // www.hyperelliptic.org / EFD / g1p / auto - shortw - projective.html  # addition-add-1998-cmo-2
    # TODO: IS THIS COMPLETE? DO WE NEED DOUBLING??
    def __add__(self, other):
        assert type(self) == type(other)
        x1, y1, z1 = self.x, self.y, self.z
        x2, y2, z2 = other.x, other.y, other.z
        y1z2 = y1 * z2
        x1z2 = x1 * z2
        z1z2 = z1 * z2
        u = y2 * z1 - y1z2
        uu = u**2
        v = x2 * z1 - x1z2
        vv = v**2
        vvv = v * vv
        R = vv * x1z2
        A = uu * z1z2 - vvv - (R + R)  # was 2*R
        x3 = v * A
        y3 = u * (R - A) - vvv * y1z2
        z3 = vvv * z1z2
        return type(self)(x3, y3, z3)

    # See https://eprint.iacr.org/2015/1060.pdf Algorithm 7 p12
    def xx__add__(self, other):
        assert type(self) == type(other)
        x1, y1, z1 = self.x, self.y, self.z
        x2, y2, z2 = other.x, other.y, other.z
        t0 = x1 * x2        # 1.  t0 ← X1 · X2
        t1 = y1 * y2        # 2.  t1 ← Y1 · Y2
        t2 = z1 * z2        # 3.  t2 ← Z1 · Z2
        t3 = x1 + y1        # 4.  t3 ← X1 + Y1
        t4 = x2 + y2        # 5.  t4 ← X2 + Y2
        t3 = t3 * t4        # 6.  t3 ← t3 · t4
        t4 = t0 + t1        # 7.  t4 ← t0 + t1
        t3 = t3 - t4        # 8.  t3 ← t3 − t4
        t4 = y1 + z1        # 9.  t4 ← Y1 + Z1
        x3 = y2 + z2        # 10. X3 ← Y2 + Z2
        t4 = t4 * x3        # 11. t4 ← t4 · X3
        x3 = t1 + t2        # 12. X3 ← t1 + t2
        t4 = t4 - x3        # 13. t4 ← t4 − X3
        x3 = x1 + z1        # 14. X3 ← X1 + Z1
        y3 = x2 + z2        # 15. Y3 ← X2 + Z2
        x3 = x3 * y3        # 16. X3 ← X3 · Y3
        y3 = t0 + t2        # 17. Y3 ← t0 + t2
        y3 = x3 - y3        # 18. Y3 ← X3 − Y3
        x3 = t0 + t0        # 19. X3 ← t0 + t0
        t0 = x3 + t0        # 20. t0 ← X3 + t0
        t2 = self.b3 * t2   # 21. t2 ← b3 · t2
        z3 = t1 + t2        # 22. Z3 ← t1 + t2
        t1 = t1 - t2        # 23. t1 ← t1 − t2
        y3 = self.b3 * y3   # 24. Y3 ← b3 · Y3
        x3 = t4 * y3        # 25. X3 ← t4 · Y3
        t2 = t3 * t1        # 26. t2 ← t3 · t1
        x3 = t2 - x3        # 27. X3 ← t2 − X3
        y3 = y3 * t0        # 28. Y3 ← Y3 · t0
        t1 = t1 * z3        # 29. t1 ← t1 · Z3
        y3 = t1 + y3        # 30. Y3 ← t1 + Y3
        t0 = t0 * t3        # 31. t0 ← t0 · t3
        z3 = z3 * t4        # 32. Z3 ← Z3 · t4
        z3 = z3 + t0        # 33. Z3 ← Z3 + t0
        return type(self)(x3, y3, z3)

    def __repr__(self):
        if self.z == type(self.z)(0): return "Neutral point"
        return f'{self.__class__.__name__} x={self.x/self.z} y={self.y/self.z}'

    def __mul__(self, other):
        assert (type(self) == Pallas and type(other) == Fq) or (type(self) == Vesta and type(other) == Fp)
        result = type(self).neutral()
        scalar = other.value
        pp = self
        while scalar > 0:
            if scalar & 0x1 != 0: result = result + pp
            scalar = scalar >> 1
            pp = pp + pp
        return result

# y**2 = x**3 + 5 over Fp
class Pallas(Curve):
    b3 = Fp(5)
    order = Fq.modulus

    @staticmethod
    def base():
        return Pallas(Fp(1), Fp(0x248b4a5cf5ed6c83ac20560f9c8711ab92e13d27d60fb1aa7f5db6c93512d546), Fp(1))

    @staticmethod
    def neutral():
        return Pallas(Fp(0), Fp(1), Fp(0))

# y**2 = x**3 + 5 over Fq
class Vesta(Curve):
    b3 = Fq(5)
    order = Fp.modulus

    @staticmethod
    def base():
        return Vesta(Fq(1), Fq(0x26bc999156dd5194ec49b1c551768ab375785e7ce00906d13e0361674fd8959f), Fq(1))

    @staticmethod
    def neutral():
        return Vesta(Fq(0), Fq(1), Fq(0))


if __name__ == "__main__":
    print("Starting quick self-test")
    from random import randint

    p1 = Pallas.neutral()
    for _i in range(100):
        p1 = p1 + Pallas.base()
    assert p1 == Pallas.base() * Fq(100)

    assert Pallas.base() * Fq(Pallas.order) == Pallas.neutral()

    v1 = Vesta.neutral()
    for _i in range(100):
        v1 = v1 + Vesta.base()
    assert v1 == Vesta.base() * Fp(100)

    for _i in range(100):
        x1 = randint(2, Pallas.order // 2)
        p1 = Pallas.base() * Fq(x1)

        x2 = randint(2, Pallas.order // 2)
        p2 = Pallas.base() * Fq(x2)

        p3 = Pallas.base() * Fq(x1 + x2)
        assert p1 + p2 == p3

    for _i in range(100):
        x1 = randint(2, Vesta.order // 2)
        v1 = Vesta.base() * Fp(x1)

        x2 = randint(2, Pallas.order // 2)
        v2 = Vesta.base() * Fp(x2)

        v3 = Vesta.base() * Fp(x1 + x2)
        assert v1 + v2 == v3

    print("Success.")
