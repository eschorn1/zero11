# This code implements the Fp and Fq fields (and helper functions) for Vesta and Pallas respectively

import secrets
from hashlib import blake2b


# Field is not meant to be used directly; it is subclassed for Fp and Fq
class __Field:
    # Constants are to be defined in specific Field subclasses
    modulus = c = s = q = None

    def __init__(self, value):
        assert type(value) is int
        self.value = value % self.modulus

    def __add__(self, other):
        assert type(self) is type(other)
        return type(self)(self.value + other.value)

    def __sub__(self, other):
        assert type(self) is type(other)
        return type(self)(self.value - other.value)

    def __mul__(self, other):
        if type(self) is not type(other): return other.__mul__(self)
        return type(self)(self.value * other.value)

    def __eq__(self, other):
        assert type(self) is type(other)
        return self.value == other.value

    def __pow__(self, other):
        assert type(other) is int
        return type(self)(pow(self.value, other, self.modulus))

    def __truediv__(self, other):
        assert type(self) is type(other)
        return self * other.inv0()

    def __neg__(self):
        return type(self)(0 - self.value)

    def __repr__(self):
        return f'{self.__class__.__name__} v={hex(self.value)}'

    def inv0(self):  # Multiplicative inverse via Fermat's little theorem (0 -> 0)
        return type(self)(pow(self.value, self.modulus - 2, self.modulus))

    def is_square(self):
        legendre_symbol = self ** ((self.modulus - 1) // 2)
        return legendre_symbol == type(self)(0) or legendre_symbol == type(self)(1)

    # Tonelli-Shanks, see https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
    #   and/or https://www.diva-portal.org/smash/get/diva2:1581080/FULLTEXT01.pdf
    def sqrt(self):
        assert self.is_square()
        m, c, i = self.s, self.c, 0
        t = pow(self.value, self.q, self.modulus)
        r = pow(self.value, (self.q + 1) // 2, self.modulus)
        while True:
            if t == 1: return type(self)(r)
            tc = t
            for i in range(m):
                if tc == 1: break
                tc = (tc * tc) % self.modulus
            b = pow(c, 2 ** (m - i - 1), self.modulus)
            m = i
            c = pow(b, 2, self.modulus)
            t = (t * c) % self.modulus
            r = (r * b) % self.modulus

    def sgn0(self):
        return self.value & 0x01  # The 'sign' of the field element

    @classmethod
    def rnd(cls):
        return cls(secrets.randbits(256 + 64))  # Oversampled to remove bias

    # Similar to https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-13.html#name-hash_to_field-implementatio
    # Code follows https://github.com/zcash/pasta_curves/blob/738fb60796d39b33f5b0a0337b8fabfcc81f98e1/src/hashtocurve.rs#L10-L77
    @classmethod
    def hash_to_field(cls, domain_prefix: bytes, message: bytes):
        curve_id = b'pallas' if cls == Fp else b'vesta'
        suffix = domain_prefix + b'-' + curve_id + b'_XMD:BLAKE2b_SSWU_RO_' + \
            bytes([22 + len(curve_id) + len(domain_prefix)])
        hasher0 = blake2b(digest_size=64, person=b'\x00' * 16)
        hasher0.update(b'\x00' * 128 + message + b'\x00\x80\x00' + suffix)
        hasher1 = blake2b(digest_size=64, person=b'\x00' * 16)
        hasher1.update(hasher0.digest() + b'\x01' + suffix)
        hasher2 = blake2b(digest_size=64, person=b'\x00' * 16)
        hasher2.update(bytes(a ^ b for (a, b) in zip(hasher0.digest(), hasher1.digest())))
        hasher2.update(b'\x02' + suffix)
        element0 = cls(int.from_bytes(hasher1.digest(), byteorder='big'))
        element1 = cls(int.from_bytes(hasher2.digest(), byteorder='big'))
        return element0, element1

    @classmethod
    def neutral(cls):
        return cls(0)


class Fp(__Field):
    s = 32  # write modulus = 2**s * q + 1 where q is odd
    q = 0x40000000000000000000000000000000224698fc094cf91b992d30ed
    modulus = 2 ** s * q + 1  # modulus = 1 mod 4, thus Tonelli-Shanks
    n = 5  # First non-square element
    c = pow(n, q, modulus)  # For sqrt()


class Fq(__Field):
    s = 32  # write modulus = 2**s * q + 1 where q is odd
    q = 0x40000000000000000000000000000000224698fc0994a8dd8c46eb21
    modulus = 2 ** s * q + 1  # modulus = 1 mod 4, thus Tonelli-Shanks
    n = 5  # First non-square element
    c = pow(n, q, modulus)  # For sqrt()


if __name__ == "__main__":
    print("Starting fields.py quick self-test")

    pp1 = Fp(3)
    pp2 = pp1 + pp1
    zp1 = pp2 / pp1
    zp2 = zp1 - Fp(1)
    assert zp2 == Fp(1)

    qq1 = Fq(3)
    qq2 = qq1 + qq1
    zq1 = qq2 / qq1
    zq2 = zq1 - Fq(1)
    assert zq2 == Fq(1)

    for _i in range(100):
        x = secrets.randbelow(Fp.modulus - 2)
        sqrt = (Fp(x) * Fp(x)).sqrt()
        assert sqrt == Fp(x) or sqrt == Fp(-x)

    for _i in range(100):
        x = secrets.randbelow(Fq.modulus - 2)
        sqrt = (Fq(x) * Fq(x)).sqrt()
        assert sqrt == Fq(x) or sqrt == Fq(-x)

    print("Success.")
