import numpy as np
import os
import pylab as pl
from core.diffract import Diffract
from core.symmetry import Symmetry
from core.atom import Atom

class Lattice(object):
    def __init__(self, pars, sym):
        self.alpha = np.deg2rad(90)
        self.beta  = np.deg2rad(90)
        self.gamma = np.deg2rad(90)
        self.sym   = sym
        self.atom     = []
        self.Symmetry = None
        self.diffract = None

        self.genSymmetry()

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

    def __repr__(self):
        s  = str([self.a,self.b,self.c])+"; "
        s += str(np.rad2deg([self.alpha,self.beta,self.gamma]))
        s += "\nSpacegroup: "+self.sym
        s += "\nAtoms:\n"
        for at in self.atom:
            s += str(at)+"\n"
        return s

    # Attach symmetry to lattice
    def genSymmetry(self):
        self.Symmetry = Symmetry(self.sym)

    # Attach diffraction to lattice
    def diffraction(self, energy, polarization=[1,0], surface=None, aziref=[0,1,0], absorb=False):
        self.diffract = Diffract(self, energy, polarization, surface, aziref, absorb)

    def genTensor(self,r,q):
        """
        ATS tensor (wrapped from symmetry class)

        :param r: position of resonant atom
        :param q: scattering vctor
        :return: ATS Tensor

        """
        if self.Symmetry == None:
            self.genSymmetry()
        return self.Symmetry.genTensor(r,q)

    def addAtom(self, element, r):
        fn = os.path.join(os.path.dirname(__file__),'dat_files'+os.sep+'atom.dat')
        f = open(fn,'r')

        Z = 0
        for line in f:
            Z += 1
            if element in line[0:2]:
                break
        f.close()
        rPos, _ = self.Symmetry.genPos(r)

        self.atom.append(Atom(element,r,rPos))

    def qMult(self,q):
        qTemp = np.dot(self.Symmetry.symOp, q)
        qEqv  = np.vstack({tuple(row) for row in qTemp})
        m = len(qEqv)

        return m, qEqv

    def genSF(self, Q):
        if(np.sum(Q)==0):
            return None
        F = np.zeros(len(self.atom),'complex128')

        for i, at in enumerate(self.atom):
            F[i] = self.chop( (np.exp(2j*np.pi*np.einsum('j,ij->i', Q, at.pos))).sum(0) )
        return F

    def chop(self,a):
        tol = 1e-12
        b=complex(a.real[abs(a.real) > tol] or 0,
                  a.imag[abs(a.imag) > tol] or 0)
        return b