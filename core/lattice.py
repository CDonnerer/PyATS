import numpy as np

from core.diffract import Diffract
from core.symmetry import Symmetry


class Lattice(object):
    def __init__(self, pars, sym):
        self.alpha = np.deg2rad(90)
        self.beta  = np.deg2rad(90)
        self.gamma = np.deg2rad(90)
        self.sym = sym
        self.pos      = None
        self.Symmetry = None
        self.diffract = None

        if len(pars)==1:
            self.a=pars[0]
            self.b=self.a
            self.c=self.a

        elif len(pars)==3:
            self.a = pars[0]
            self.b = pars[1]
            self.c = pars[2]

        else:
            self.a = pars[0]
            self.b = pars[1]
            self.c = pars[2]
            self.alpha = np.deg2rad(pars[3])
            self.beta  = np.deg2rad(pars[4])
            self.gamma = np.deg2rad(pars[5])

    # Attach symmetry to lattice
    def genSymmetry(self):
        self.Symmetry = Symmetry(self.sym)

    # Attach diffraction to lattice
    def diffraction(self, energy, polarization, surface, aziref):
        self.diffract = Diffract(self, energy, polarization, surface, aziref)

    # ATS tensor (wrapped from symmetry)
    def genTensor(self,r,q):
        if self.Symmetry == None:
            self.genSymmetry()
        return self.Symmetry.genTensor(r,q)

    # XRMS tensor
    def xrmsTensor(self, M):
        Fm = 1j * np.array([[0, M[2], -M[1]], [-M[2], 0, M[0]], [M[1], -M[0], 0]])
        return Fm

    def addAtom(self, r, element):
        f = open('dat_files/atom.dat')
        Z = 0
        for line in f:
            Z += 1
            if element in line[0:2]:
                #print(Z," ",line, " at ", r)
                break
        f.close()
        rPos, _ = self.Symmetry.genPos(r)
        self.pos = rPos

    def qMult(self,q):
        qTemp = np.dot(self.lat.Symmetry.symOp, q)
        qEqv = np.vstack({tuple(row) for row in qTemp})
        m = len(qEqv)

        return m, qEqv

    def genSF(self, Q):
        F = complex(0,0)

        for atom in self.pos:
            F += np.exp(2j * np.pi * np.dot(Q, atom.T))

        return self.chop(F)

    def chop(self,a):
        tol = 1e-8
        b=complex(a.real[abs(a.real) > tol] or 0, a.imag[abs(a.imag) > tol] or 0)
        return b