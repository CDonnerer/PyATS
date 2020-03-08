import os
import numpy as np

class Symmetry(object):
    def __init__(self, sym):
        self.sym = sym
        self.symOp, self.symTr = self.generate_operators(sym)

    def generate_operators(self,sym):
        """
        Generate spacegroup symmetry operations

        :param sym:
        :return: symOp, symTr:
        """

        tol = 1e-5
        gener, trans, symStr = self.generate_symmetry(sym)

        nGen = len(gener[1, 1, :])
        symOp = []
        symTr = []
        symOp.append(np.identity(3))
        symTr.append(np.zeros(3))
        P = np.zeros(nGen, dtype=np.int)

        for i in range(0, nGen):
            R0 = gener[:, :, i]
            T0 = trans[:, i]
            P[i] = self.symmetry_order(R0, T0)

            R = np.identity(3)
            T = np.zeros(3)
            for j in range(0, P[i] - 1):
                R = np.dot(R0, R)
                T = R0.dot(T) + T0
                nSym = len(symTr)

                for k in range(0, nSym):
                    RS = np.dot(R, symOp[k])
                    TS = (R.dot(symTr[k]) + T) % 1

                    idxR = np.sum(np.sum(abs(symOp - RS), axis=1), axis=1) > tol
                    idxT = np.sum(abs(symTr - TS), axis=1) > tol

                    if all(idxT | idxR):
                        symTr.append(TS)
                        symOp.append(RS)
        return symOp, symTr

    def generate_symmetry(self,sym):
        """
        Read in spacegoup generators from file

        :param sym:
        :return: symOp, symTr, symStr:
        """
        fn = os.path.join(os.path.dirname(__file__),'dat_files'+os.sep+'symmetry.dat')
        f = open(fn,'r')

        for line in f:
            if sym in line:
                symStr = line[19:]
                break
        f.close()

        vNew = np.zeros(3)
        symOp = np.zeros(shape=(3, 3, 20))
        symTr = np.zeros(shape=(3, 20))
        nNew = 0
        nOp = 0
        nSign = 1

        j = 0
        while (j < len(symStr)):
            if symStr[j] == ',':
                symOp[nNew, :, nOp] = vNew
                vNew *= 0
                nSign = 1
                nNew = nNew % 2 + 1
            elif symStr[j] == ';':
                symOp[nNew, :, nOp] = vNew
                vNew *= 0
                nSign = 1
                nNew = 0
                nOp += 1
            elif symStr[j] == 'x':
                vNew[0] = nSign
            elif symStr[j] == 'y':
                vNew[1] = nSign
            elif symStr[j] == 'z':
                vNew[2] = nSign
            elif symStr[j] == '-':
                nSign = -1
            elif symStr[j] == '+':
                nSign = 1
            elif symStr[j] == '1' or symStr[j] == '2' or symStr[j] == '3':
                symTr[nNew, nOp] = float(symStr[j]) / float(symStr[j + 2])
                j += 2
            j += 1

        symOp[nNew, :, nOp] = vNew
        symOp = symOp[:, :, 0:nOp+1]
        symTr = symTr[:, 0:nOp+1]
        return symOp, symTr, symStr

    def symmetry_order(self,R, T):
        """
        Determine order of symmetry operator

        :param R:
        :param T:
        :return:
        """
        tol = 1e-5

        N  = 1
        RN = R
        TN = T

        while (np.linalg.norm(RN - np.identity(3)) > tol or np.linalg.norm(TN) > tol) and N < 10:
            RN = np.dot(R, RN)
            TN = (R.dot(TN) + T) % 1
            N += 1

        return N

    def generate_positions(self, r):
        """
        Generate symmetry-equivalent positions

        :param r:
        :return: rPos, isMoved:
        """
        rTemp = (np.dot(self.symOp, r) + self.symTr) % 1
        rPos = np.vstack({tuple(row) for row in rTemp})

        isMoved = np.sum(rTemp - r, axis=1) > 0

        return rPos, isMoved

    def point_symmetry(self, r):
        """
        Generate point group symmetry operations

        :param r:
        :return: pointOp:
        """
        _, isMoved = self.generate_positions(r)

        opArray = np.asarray(self.symOp)
        pointOp = opArray[np.logical_not(isMoved)]

        return pointOp

    def generate_tensor(self, r, q):
        """
        Generate ATS scattering tensor

        :param r: site of absorbing atom
        :param q: momentum transfer
        :return: 3x3 scattering tensor
        """

        anisoT = np.ones(shape=(3, 3)) - np.identity(3)

        b  = np.exp(2j * np.pi * np.einsum('i,ji->j', q, (self.symTr + np.dot(self.symOp, r))))
        a  = np.einsum('mij,jk,mlk->mil', self.symOp, anisoT, self.symOp)
        Fk = np.einsum('ijk,i->jk', a, b)

        self.chop(Fk)
        norm = Fk.sum() / 2

        # Check if tensor is non-zero
        with np.warnings.catch_warnings():
            np.warnings.filterwarnings('error')
            try:
                return Fk / norm
            except:
                print('Reflection Q = (', q[0], q[1], q[2],') is not ATS allowed!')
                return None

    def chop(self,a):
        """
        Helper function to get rid of numerical errors

        :param a:
        :return:
        """
        tol = 1e-8
        a.real[abs(a.real) < tol] = 0
        a.imag[abs(a.imag) < tol] = 0
        return a