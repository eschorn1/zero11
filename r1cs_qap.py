from fields import Fp
from poly import Poly


class R1csQap:

    def __init__(self, input_names):
        assert isinstance(input_names, list) and (len(input_names) > 0)
        self.input_names = input_names
        self.gates = dict()
        self.gates['left'] = []; self.gates['right'] = []; self.gates['out'] = []
        self.samples = dict()
        self.samples['left'] = []; self.samples['right'] = []; self.samples['out'] = []
        self.polys = dict()
        self.polys['left'] = []; self.polys['right'] = []; self.polys['out'] = []
        self.t = None; self.h = None; self.z = None

    def append_gate(self, left_tuple, right_tuple, output_tuple):  # [(name, int value)]
        left_list = [Fp(0)] * len(self.input_names)
        right_list = [Fp(0)] * len(self.input_names)
        output_list = [Fp(0)] * len(self.input_names)
        for item in left_tuple: left_list[self.input_names.index(item[0])] = Fp(item[1])
        for item in right_tuple: right_list[self.input_names.index(item[0])] = Fp(item[1])
        for item in output_tuple: output_list[self.input_names.index(item[0])] = Fp(item[1])
        self.gates['left'].append(left_list); self.gates['right'].append(right_list); self.gates['out'].append(output_list)

    def transpose(self):
        for matrix in ['left', 'right', 'out']:
            rows = len(self.gates['left'])
            columns = len(self.gates['left'][0])
            result = []
            for j in range(columns):
                row = []
                for i in range(rows):
                    row.append(self.gates[matrix][i][j])
                result.append(row)
            self.samples[matrix] = result

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

    def gen_polys(self):
        for group in ['left', 'right', 'out']:
            for samples in self.samples[group]:
                res_poly = self.poly_interp(samples)
                self.polys[group].append(res_poly)

    def gen_t(self, soln):
        dots = dict()
        for group in ['left', 'right', 'out']:
            xxx = Poly([Fp(0)])
            for index in range(len(soln)):
                xxx = xxx + self.polys[group][index] * Poly([soln[index]])
            dots[group] = xxx
        self.t = (dots['left'] * dots['right']) - dots['out']

    def gen_z(self):
        z = Poly([Fp(1)])
        for index in range(1, 5):
            z = z * Poly([Fp(-index), Fp(1)])
        self.z = z


if __name__ == '__main__':

    gates = R1csQap(['one', 'x', 'out', 'sym_1', 'y', 'sym_2'])
    gates.append_gate([('x', 1)], [('x', 1)], [('sym_1', 1)])
    gates.append_gate([('sym_1', 1)], [('x', 1)], [('y', 1)])
    gates.append_gate([('x', 1), ('y', 1)], [('one', 1)], [('sym_2', 1)])
    gates.append_gate([('one', 5), ('sym_2', 1)], [('one', 1)], [('out', 1)])    # PRINT

    for gate in gates.gates['left']: print(f'A gate {gate}')
    for gate in gates.gates['right']: print(f'B gate {gate}')
    for gate in gates.gates['out']: print(f'C gate {gate}')

    gates.transpose()

    for sample in gates.samples['left']: print(f'A samp {sample}')
    for sample in gates.samples['right']: print(f'B samp {sample}')
    for sample in gates.samples['out']: print(f'C samp {sample}')

    gates.gen_polys()

    for poly in gates.polys['left']: print(f'A poly {poly}')
    for poly in gates.polys['right']: print(f'B poly {poly}')
    for poly in gates.polys['out']: print(f'C poly {poly}')

    soln1 = [Fp(1), Fp(3), Fp(35), Fp(9), Fp(27), Fp(30)]
    gates.gen_t(soln1)

    gates.gen_z()

    print(" h is ", gates.t / gates. z)
