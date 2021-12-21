from fields import Fp
from poly import Poly

class Gates:
    def __init__(self, input_names):
        self.inputs = input_names
        self.a = []
        self.b = []
        self.c = []

    def append(self, left, right, output):  # [(name, value)]
        l_list = [Fp(0)] * len(self.inputs)
        r_list = [Fp(0)] * len(self.inputs)
        o_list = [Fp(0)] * len(self.inputs)
        for item in left: l_list[self.inputs.index(item[0])] = Fp(item[1])
        for item in right: r_list[self.inputs.index(item[0])] = Fp(item[1])
        for item in output: o_list[self.inputs.index(item[0])] = Fp(item[1])
        self.a.append(l_list); self.b.append(r_list); self.c.append(o_list)

    def transpose(self):
        result = []
        for matrix in [self.a, self.b, self.c]:
            rows = len(self.a)
            columns = len(self.a[0])
            res = []
            for j in range(columns):
                row = []
                for i in range(rows):
                    row.append(matrix[i][j])
                res.append(row)
            result.append(res)
        return result

    def pprint(self):
        print("A")
        for gate in self.a: print(gate)
        print("B")
        for gate in self.b: print(gate)
        print("C")
        for gate in self.c: print(gate)

    @staticmethod
    def poly_interp(samples):
        result = Poly([Fp(0)])
        for (sample_index, value) in enumerate(samples):
            terms = Poly([Fp(1)])
            if value == Fp(0): continue
            zero_indices = [x for x in range(len(samples))]
            zero_indices.remove(sample_index)
            for zero_index in zero_indices:
                terms = terms * Poly([Fp(-(zero_index + 1)), Fp(1)])  # X - index (starting at one)
            scale = terms.eval(Fp(sample_index + 1))
            addd = Poly([value]) * (terms / Poly([scale]))[0]
            result = result + addd
        return result


if __name__ == '__main__':

    # Build R1CS; TODO: maybe can use a multidim dict with default = 0 since sparse
    gates = Gates(['one', 'x', 'out', 'sym_1', 'y', 'sym_2'])
    gates.append([('x', 1)], [('x', 1)], [('sym_1', 1)])
    gates.append([('sym_1', 1)], [('x', 1)], [('y', 1)])
    gates.append([('x', 1), ('y', 1)], [('one', 1)], [('sym_2', 1)])
    gates.append([('one', 5), ('sym_2', 1)], [('one', 1)], [('out', 1)])
    # PRINT
    gates.pprint()

    # Confirm each a * b - c == 0
    soln = [Fp(1), Fp(3), Fp(35), Fp(9), Fp(27), Fp(30)]
    for index in range(len(gates.a)):
        a = Poly.dot(soln, gates.a[index])
        b = Poly.dot(soln, gates.b[index])
        c = Poly.dot(soln, gates.c[index])
        print(a * b - c)
        assert(a * b - c == Fp(0))

    zz = gates.transpose()

    polys = []
    for group in zz:
        poly_group = []
        for entry in group:
            res_poly = Gates.poly_interp(entry)
            poly_group.append(res_poly)
            for val in range(1, 5):
                print("INTERP", val, res_poly.eval(Fp(val)), res_poly)
        polys.append(poly_group)

    print("EVAL A POLYS")
    for index, poli in enumerate(polys[1]):
        print(4, poli.eval(Fp(4)))


    dots = []
    for group in polys:
        xxx = Poly([Fp(0)])
        for index in range(len(soln)):
            xxx = xxx + group[index] * Poly([soln[index]])
#        print("xxx", xxx)
        dots.append(xxx)

    print("dots", dots)

    t = (dots[0] * dots[1]) - dots[2]
    print("ttt", t)

    for index in range(1, 5):
        print("t eval", index, t.eval(Fp(index)))


    z = Poly([Fp(1)])
    for index in range(1, 5):
        z = z * Poly([Fp(-index), Fp(1)])

    h = t / z

    print("hhh", h)