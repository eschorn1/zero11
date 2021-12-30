import secrets

from fields import Fp, Fq


class Curve:

    # Constants are to be defined in specific Field subclasses
    A = B = Z = theta = iso = a = b3 = None  # Code will fail if Field is directly utilized

    def __init__(self, x, y, z):
        assert type(x) == Fp or type(x) == Fq
        assert type(x) == type(y) and type(y) == type(z)
        self.x, self.y, self.z = x, y, z

    def __eq__(self, other):
        assert type(self) == type(other)
        if self.z == (type(self.z))(0): return other.z == (type(other.z))(0)
        if other.z == (type(other.z))(0): return False
        return self.x/self.z == other.x/other.z and self.y/self.z == other.y/other.z

    # See https://eprint.iacr.org/2015/1060.pdf page 8  # TODO: Check b3 = 3*b
    # Algorithm 1: Complete, projective point addition for arbitrary prime order short Weierstrass curves E/Fq : y^2 = x^3 + ax + b.
    def __add__(self, other, a=None, b3=None):
        assert type(self) == type(other)
        # if a is None and b3 is None: assert self.is_on_curve(self)
        x1, y1, z1 = self.x, self.y, self.z
        x2, y2, z2 = other.x, other.y, other.z
        if a is None: a = self.a
        if b3 is None: b3 = self.b3 * type(self.x)(3)
        assert self.is_on_my_curve(a, b3 / type(self.x)(3))
        t0 = x1 * x2     # 1.  t0 ← X1 · X2
        t1 = y1 * y2     # 2.  t1 ← Y1 · Y2
        t2 = z1 * z2     # 3.  t2 ← Z1 · Z2
        t3 = x1 + y1     # 4.  t3 ← X1 + Y1
        t4 = x2 + y2     # 5.  t4 ← X2 + Y2
        t3 = t3 * t4     # 6.  t3 ← t3 · t4
        t4 = t0 + t1     # 7.  t4 ← t0 + t1
        t3 = t3 - t4     # 8.  t3 ← t3 − t4
        t4 = x1 + z1     # 9.  t4 ← X1 + Z1
        t5 = x2 + z2     # 10. t5 ← X2 + Z2
        t4 = t4 * t5     # 11. t4 ← t4 · t5
        t5 = t0 + t2     # 12. t5 ← t0 + t2
        t4 = t4 - t5     # 13. t4 ← t4 − t5
        t5 = y1 + z1     # 14. t5 ← Y1 + Z1
        x3 = y2 + z2     # 15. X3 ← Y2 + Z2
        t5 = t5 * x3     # 16. t5 ← t5 · X3
        x3 = t1 + t2     # 17. X3 ← t1 + t2
        t5 = t5 - x3     # 18. t5 ← t5 − X3
        z3 = a * t4      # 19. Z3 ← a · t4
        x3 = b3 * t2     # 20. X3 ← b3 · t2
        z3 = x3 + z3     # 21. Z3 ← X3 + Z3
        x3 = t1 - z3     # 22. X3 ← t1 − Z3
        z3 = t1 + z3     # 23. Z3 ← t1 + Z3
        y3 = x3 * z3     # 24. Y3 ← X3 · Z3
        t1 = t0 + t0     # 25. t1 ← t0 + t0
        t1 = t1 + t0     # 26. t1 ← t1 + t0
        t2 = a * t2      # 27. t2 ← a · t2
        t4 = b3 * t4     # 28. t4 ← b3 · t4
        t1 = t1 + t2     # 29. t1 ← t1 + t2
        t2 = t0 - t2     # 30. t2 ← t0 − t2
        t2 = a * t2      # 31. t2 ← a · t2
        t4 = t4 + t2     # 32. t4 ← t4 + t2
        t0 = t1 * t4     # 33. t0 ← t1 · t4
        y3 = y3 + t0     # 34. Y3 ← Y3 + t0
        t0 = t5 * t4     # 35. t0 ← t5 · t4
        x3 = t3 * x3     # 36. X3 ← t3 · X3
        x3 = x3 - t0     # 37. X3 ← X3 − t0
        t0 = t3 * t1     # 38. t0 ← t3 · t1
        z3 = t5 * z3     # 39. Z3 ← t5 · Z3
        z3 = z3 + t0     # 40. Z3 ← Z3 + t0
        result = type(self)(x3, y3, z3)
        # assert self.is_on_curve(result)
        return result

    def __repr__(self):
        if self.z == (type(self.z))(0): return "Neutral point"
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

    @classmethod
    def hash_to_curve(cls, message):
        e1, e2 = type(cls.b3).hash_to_field(cls.__name__.lower().encode('utf-8'), b'z.cash:test', message)
        (x1, y1) = cls.map_to_curve_simple_swu(e1)
        (x2, y2) = cls.map_to_curve_simple_swu(e2)
        q0 = cls(x1, y1, type(cls.b3)(1))
        q1 = cls(x2, y2, type(cls.b3)(1))
        # r = q0 + q1
        r = q0.__add__(q1, a=cls.A, b3=cls.B * Fp(3))
        assert r.is_on_my_curve(a=cls.A, b=cls.B * Fp(1))
        z = cls.iso_map(r)
        assert z.is_on_my_curve(a=cls.a, b=cls.b3 * Fp(1))
        return z

    @classmethod
    def map_to_curve_simple_swu(cls, u):  # TODO: This should probably return type iso_<curve>
        tv1 = type(u)(1) / (cls.Z ** 2 * u ** 4 + cls.Z * u ** 2)  # 1. tv1 = inv0(Z^2 * u^4 + Z * u^2)
        x1 = ((type(u)(0) - cls.B) / cls.A) * (type(u)(1) + tv1)  # 2.  x1 = (-B / A) * (1 + tv1)
        if tv1 == type(u)(0): x1 = cls.B / (cls.Z * cls.A)  # 3.  If tv1 == 0, set x1 = B / (Z * A)
        gx1 = x1 ** 3 + cls.A * x1 + cls.B  # 4. gx1 = x1^3 + A * x1 + B
        x2 = cls.Z * u ** 2 * x1  # 5.  x2 = Z * u^2 * x1
        gx2 = x2 ** 3 + cls.A * x2 + cls.B  # 6. gx2 = x2^3 + A * x2 + B
        if gx1.is_square():  # 7.  If is_square(gx1), set x = x1 and y = sqrt(gx1)
            x = x1
            y = gx1.sqrt()
        else:  # 8.  Else set x = x2 and y = sqrt(gx2)
            x = x2
            y = gx2.sqrt()
        if u.sgn0 != y.sgn0: y = type(u)(0) - y  # 9.  If sgn0(u) != sgn0(y), set y = -y
        return x, y  # 10. return (x, y)

    def is_on_my_curve(self, a, b):
        return self.z * self.y ** 2 == self.x ** 3 + a * self.x * self.z ** 2 + b * self.z ** 3

    @classmethod
    def iso_map(cls, pt):
        x = pt.x / pt.z
        y = pt.y / pt.z
        x_num = cls.iso[0] * x ** 3 + cls.iso[1] * x ** 2 + cls.iso[2] * x + cls.iso[3]
        x_den = x ** 2 + cls.iso[4] * x + cls.iso[5]
        y_num = cls.iso[6] * x ** 3 + cls.iso[7] * x ** 2 + cls.iso[8] * x + cls.iso[9]
        y_den = x ** 3 + cls.iso[10] * x ** 2 + cls.iso[11] * x + cls.iso[12]
        return type(pt)(x_num / x_den, (y * y_num) / y_den, Fp(1))

    @classmethod
    def rnd(cls):
        token_g = secrets.token_bytes(16)
        element = Pallas.hash_to_curve(token_g)
        return element

    @classmethod
    def neutral(cls):
        pass


# y**2 = x**3 + 5 over Fp
class Pallas(Curve):
    a = Fp(0)
    b3 = Fp(5)
    order = Fq.modulus
    #     "iso-pallas",
    A = Fp(0x18354a2eb0ea8c9c49be2d7258370742b74134581a27a59f92bb4b0b657a014b)
    B = Fp(1265)
    Z = Fp(-13)
    theta = Fp(0x0f7bdb65814179b44647aef782d5cdc851f64fc4dc888857ca330bcc09ac318e)
    iso = [Fp(0x0e38e38e38e38e38e38e38e38e38e38e4081775473d8375b775f6034aaaaaaab),
           Fp(0x3509afd51872d88e267c7ffa51cf412a0f93b82ee4b994958cf863b02814fb76),
           Fp(0x17329b9ec525375398c7d7ac3d98fd13380af066cfeb6d690eb64faef37ea4f7),
           Fp(0x1c71c71c71c71c71c71c71c71c71c71c8102eea8e7b06eb6eebec06955555580),
           Fp(0x1d572e7ddc099cff5a607fcce0494a799c434ac1c96b6980c47f2ab668bcd71f),
           Fp(0x325669becaecd5d11d13bf2a7f22b105b4abf9fb9a1fc81c2aa3af1eae5b6604),
           Fp(0x1a12f684bda12f684bda12f684bda12f7642b01ad461bad25ad985b5e38e38e4),
           Fp(0x1a84d7ea8c396c47133e3ffd28e7a09507c9dc17725cca4ac67c31d8140a7dbb),
           Fp(0x3fb98ff0d2ddcadd303216cce1db9ff11765e924f745937802e2be87d225b234),
           Fp(0x025ed097b425ed097b425ed097b425ed0ac03e8e134eb3e493e53ab371c71c4f),
           Fp(0x0c02c5bcca0e6b7f0790bfb3506defb65941a3a4a97aa1b35a28279b1d1b42ae),
           Fp(0x17033d3c60c68173573b3d7f7d681310d976bbfabbc5661d4d90ab820b12320a),
           Fp(0x40000000000000000000000000000000224698fc094cf91b992d30ecfffffde5)]

    @staticmethod
    def base():
        return Pallas(Fp(1), Fp(0x248b4a5cf5ed6c83ac20560f9c8711ab92e13d27d60fb1aa7f5db6c93512d546), Fp(1))

    @staticmethod
    def neutral():
        return Pallas(Fp(0), Fp(1), Fp(0))

    neu = (Fp(0), Fp(1), Fp(0))


# y**2 = x**3 + 5 over Fq
class Vesta(Curve):
    a = Fq(0)
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
        xx1 = randint(2, Pallas.order // 2)
        p1 = Pallas.base() * Fq(xx1)

        xx2 = randint(2, Pallas.order // 2)
        p2 = Pallas.base() * Fq(xx2)

        p3 = Pallas.base() * Fq(xx1 + xx2)
        assert p1 + p2 == p3

    for _i in range(100):
        xx1 = randint(2, Vesta.order // 2)
        v1 = Vesta.base() * Fp(xx1)

        x2 = randint(2, Pallas.order // 2)
        v2 = Vesta.base() * Fp(xx2)

        v3 = Vesta.base() * Fp(xx1 + xx2)
        assert v1 + v2 == v3

    print("Testing hash-to-curve")
    tmp1, tmp2 = Fp.hash_to_field(b'pallas', b'z.cash:test', b'Trans rights now!')
    (act_x1, act_y1) = Pallas.map_to_curve_simple_swu(tmp1)
    assert act_x1 == Fp(0x05c3482fe40155e152fdc0be06c4766b67a2b3d8d9bb64ee6137382879dc2160)
    assert act_y1 == Fp(0x3825fb730c259375175ff31b94dc36dcf031b13f3116bda725f1c98717739f1f)
    (act_x2, act_y2) = Pallas.map_to_curve_simple_swu(tmp2)
    assert act_x2 == Fp(0x2c6e5aa1a88cd76c8a9d436438d2993244bf7704e4f322a86d0890bd6cee28ab)
    assert act_y2 == Fp(0x0b20c46efea44d15e4828808c86a72789d54328635ba4274d8e9b48d9654f65b)

    qq0 = Pallas(act_x1, act_y1, Fp(1))
    qq1 = Pallas(act_x2, act_y2, Fp(1))
    # r = q0 + q1
    rr = qq0.__add__(qq1, a=Pallas.A, b3=Pallas.B * Fp(3))

    assert rr == Pallas(Fp(0x3da8497a87f06e28b7983f044f5f93575daf4806e0735700ebd79184070bb58e),
                        Fp(0x271b2b52a5e759ee28a21db1a520739b1f53a1960433d08593ec10e225dec8c0),
                        Fp(1))

    assert rr.is_on_my_curve(a=rr.A, b=rr.B * Fp(1))

    z11 = Pallas.iso_map(rr)
    assert z11 == Pallas(Fp(0x1818cda31ffdc8c3ff23df3d88c26f952340257d0f187a0236695c9b640b6bd3),
                         Fp(0x1e20888510123752166a0306332e126289f6f9a2774160395f2f1efc9b1280c),
                         Fp(1))

    zz = Pallas.hash_to_curve(b'Trans rights now!')
    assert zz == z11

    print("Success.")
