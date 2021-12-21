from fields import Fp, Fq


class Poly:
    def __init__(self, coeffs):
        assert isinstance(coeffs, list) and (len(coeffs) == 0 or (isinstance(coeffs[0], Fp) or isinstance(coeffs[0], Fq)))
        self.coeffs = self.trim_leading_0s(coeffs)  # list of Fp or Fq or empty

    def __repr__(self):
        coeffs = [x.__repr__() for x in self.coeffs]
        return f'{self.__class__.__name__} coeffs={coeffs}'

    def __eq__(self, other):
        assert type(self) == type(other)
        if len(self.coeffs) != len(other.coeffs): return False
        return all(map(lambda x, y: x == y, self.coeffs, other.coeffs))

    def __add__(self, other):
        assert type(self) == type(other)
        if len(self.coeffs) == 0: return other
        if len(other.coeffs) == 0: return self
        result = []
        length = max(len(self.coeffs), len(other.coeffs))
        for index in range(length):
            if index >= len(self.coeffs):
                result.append(other.coeffs[index])
            elif index >= len(other.coeffs):
                result.append(self.coeffs[index])
            else:
                result.append(self.coeffs[index] + other.coeffs[index])
        return Poly(self.trim_leading_0s(result))

    def __sub__(self, other):
        assert type(self) == type(other)
        if len(other.coeffs) == 0: return self
        if len(self.coeffs) == 0: self.coeffs = [type(other.coeffs[0])(0)]
        result = []
        length = max(len(self.coeffs), len(other.coeffs))
        for index in range(length):
            if index >= len(self.coeffs):
                result.append(type(self.coeffs[0])(0) - other.coeffs[index])
            elif index >= len(other.coeffs):
                result.append(self.coeffs[index])
            else:
                result.append(self.coeffs[index] - other.coeffs[index])
        return Poly(self.trim_leading_0s(result))

    def __mul__(self, other):
        assert type(self) == type(other)
        if len(self.coeffs) == 0 or len(other.coeffs) == 0: return Poly([])
        result = [type(self.coeffs[0])(0)] * (len(self.coeffs) * len(other.coeffs))
        for i in range(len(self.coeffs)):
            for j in range(len(other.coeffs)):
                result[i+j] = result[i+j] + self.coeffs[i] * other.coeffs[j]
        return Poly(self.trim_leading_0s(result))

    def __truediv__(self, other):
        assert type(self) == type(other) and len(other.coeffs) != 0
        q = Poly([])
        r = self
        while len(r.coeffs) != 0 and len(r.coeffs) >= len(other.coeffs):
            t = r.coeffs[-1] / other.coeffs[-1]  # Leading scale factor
            trailing = [type(self.coeffs[0])(0)] * (len(r.coeffs) - len(other.coeffs))
            q = q + Poly(trailing + [t])
            r = r - Poly(trailing + [t]) * other
        return q, r  # at some point, require r = []

    # Note, this can/will evolve to handle scalar * point
    # https://en.wikipedia.org/wiki/Horner%27s_method
    def eval(self, x):
        assert isinstance(x, Fp) or isinstance(x, Fq)
        result = self.coeffs[-1]
        for index in range(len(self.coeffs)-1, 0, -1):
            result = result * x
            result = result + self.coeffs[index-1]
        return result

    # This probably wants to live elsewhere
    @staticmethod
    def dot(left, right):
        # assert type(left) == type(right) and len(left.coeffs) == len(right.coeffs)
        xx = map(lambda x, y: x * y, left, right)
        return sum(xx, start=left[0].neutral())

    @staticmethod
    def trim_leading_0s(coeffs):
        assert isinstance(coeffs, list) and (len(coeffs) == 0 or (isinstance(coeffs[0], Fp) or isinstance(coeffs[0], Fq)))
        if not coeffs: return coeffs
        while coeffs[-1].value == 0:
            coeffs.pop()
            if not coeffs: return coeffs  # Note: zero poly is empty poly
        return coeffs


if __name__ == "__main__":
    print("Starting quick self-test")
    from random import randint
    from curves import Pallas

    for _i in range(100):
        p1 = Poly([Fp(randint(0, Fp.modulus - 1)) for i in range(20)])
        p2 = Poly([Fp(randint(0, 20)) for i in range(40)])  # aiming for some zeros
        pmul = p1 * p2
        assert (pmul / p1)[0] == p2
        assert (pmul / p2)[0] == p1

    for _i in range(100):
        q1 = Poly([Fq(randint(0, Fq.modulus - 1)) for i in range(20)])
        q2 = Poly([Fq(randint(0, 20)) for i in range(40)])  # aiming for some zeros
        qmul = q1 * q2
        assert (qmul / q1)[0] == q2
        assert (qmul / q2)[0] == q1

    # f(x) = 1 + 2x + 3x^2 + 4x^3; f(5) = 586
    p1 = Poly([Fp(x) for x in [1, 2, 3, 4]])
    assert p1.eval(Fp(5)) == Fp(586)

    # (1 + 2 + 3 + 4) dot [100, 100, 100, 100] = 10 * 100
    assert Poly.dot(p1.coeffs, [Fp(100), Fp(100), Fp(100), Fp(100)]) == Fp(1000)

    # (1, 2, 3, 4) dot (1, 2, 3, 4) = 1+4+9+16 = 30
    right = [Fq(1), Fq(2), Fq(3), Fq(4)]
    left = [Pallas.base() * x for x in [Fq(1), Fq(2), Fq(3), Fq(4)]]
    res = Poly.dot(left, right)
    assert res == Pallas.base() * Fq(30)

    # TODO: want to evaluate a0*x^0 + a1*x^1 + a2*x^2 + a3*x^3
    #       where a0 are in Fq and x is in {Fq or Pallas point} ?????? IS RIGHT ????? READ PAPER
    # TODO: --------------> Better review paper and see what is needed; maybe it is only dot product!!!!!
    # TODO: NNNNNNOOOOOOOOOOOOOOOOo, more likely dot must handle fp x fp, fp x g
    pp = Pallas.base()


    print("Success.")

