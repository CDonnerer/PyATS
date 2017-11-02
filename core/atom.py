import os
import io
import numpy as np

class Atom(object):
    def __init__(self, element, wpos, rpos):
        """
        :param element: Periodic symbol
        :param wpos:    Wyckoff position
        :param rpos:    Symmetry-equivalent positions
        """
        self.element = element
        self.wpos = wpos
        self.pos  = rpos
        self.a, self.b, self.c = self.read_formfactor()

    def __repr__(self):
        n = len(self.pos)
        s = str(n) + " " + self.element + " at r = " + str(self.wpos)
        return s

    def read_formfactor(self):
        """
        Read-in atomic form factor parameters from file

        :return:
        """
        fn = os.path.join(os.path.dirname(__file__), 'dat_files' + os.sep + 'formfactor.dat')

        # TODO: check file. works for Python 2,3 c.f. symmetry
        f = io.open(fn, mode="r", encoding="utf-8")
        #f = open(fn, 'r')

        for line in f:
            if self.element in line:
                ff_str = line[0:]
                break
        f.close()

        a = np.zeros(4)
        b = np.zeros(4)
        tokens = ff_str.split()

        j = 1
        for i in range(0, 4):
            a[i] = tokens[j]
            b[i] = tokens[j + 1]
            j += 2
        c = float(tokens[-1])

        return a, b, c

    def formfactor(self, sinThetaOverLam):
        ff = (self.a * np.exp(-self.b * (sinThetaOverLam ** 2))).sum() + self.c
        return ff