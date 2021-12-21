class Field:

    # Constants are to be defined in specific Field subclasses
    modulus = c = s = q = None  # Code will fail if Field is directly utilized

    def __init__(self, value):
        assert type(value) == int
        self.value = value % self.modulus

    def __repr__(self):
        return f'{self.__class__.__name__} v={hex(self.value)}'

    def __eq__(self, other):
        assert type(self) == type(other)
        return self.value == other.value

    def __add__(self, other):
        assert type(self) == type(other)
        return type(self)(self.value + other.value)

    def __sub__(self, other):
        assert type(self) == type(other)
        return type(self)(self.value - other.value)

    def __mul__(self, other):
        assert type(self) == type(other)
        return type(self)(self.value * other.value)

    def __pow__(self, other):
        assert isinstance(other, int)
        return type(self)(pow(self.value, other, self.modulus))

    def __truediv__(self, other):
        assert type(self) == type(other)
        inv = pow(other.value, self.modulus - 2, self.modulus)
        return type(self)(self.value * inv)

    def is_square(self):
        val = self**((self.modulus-1)//2)
        return val == type(self)(0) or val == type(self)(1)

    # Tonelli-Shanks, see https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm and/or https://www.diva-portal.org/smash/get/diva2:1581080/FULLTEXT01.pdf
    def sqrt(self):
        m, c = self.s, self.c
        t = pow(self.value, self.q, self.modulus)
        r = pow(self.value, (self.q + 1) // 2, self.modulus)
        while True:
            if t == 1: return type(self)(r)
            tc = t
            for i in range(m):
                if tc == 1: break
                tc = (tc*tc) % self.modulus
            b = pow(c, 2**(m-i-1), self.modulus)
            m = i
            c = pow(b, 2, self.modulus)
            t = (t*c) % self.modulus
            r = (r*b) % self.modulus

    def sgn0(self):
        return self.value & 0x01


class Fp(Field):
    s = 32  # write modulus = 2**s * q + 1 where q is odd
    q = 0x40000000000000000000000000000000224698fc094cf91b992d30ed
    modulus = 2**s * q + 1  # modulus = 1 mod 4, thus Tonelli-Shanks
    n = 5  # First non-square
    c = pow(n, q, modulus)  # For sqrt()

    @staticmethod
    def neutral():
        return Fp(0)


class Fq(Field):
    s = 32  # write modulus = 2**s * q + 1 where q is odd
    q = 0x40000000000000000000000000000000224698fc0994a8dd8c46eb21
    modulus = 2**s * q + 1  # modulus = 1 mod 4, thus Tonelli-Shanks
    n = 5  # First non-square
    c = pow(n, q, modulus)  # For sqrt()

    @staticmethod
    def neutral():
        return Fq(0)


if __name__ == "__main__":
    print("Starting quick self-test")
    from random import randint

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

    for _i in range(1000):
        x = randint(2, Fp.modulus-2)
        sqrt = (Fp(x) * Fp(x)).sqrt()
        assert sqrt == Fp(x) or sqrt == Fp(-x)

    for _i in range(1000):
        x = randint(2, Fq.modulus-2)
        sqrt = (Fq(x) * Fq(x)).sqrt()
        assert sqrt == Fq(x) or sqrt == Fq(-x)

    print("Success.")
