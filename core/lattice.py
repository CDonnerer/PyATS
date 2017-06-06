import numpy as np
from core.diffract import Diffract
from core.symmetry import Symmetry
from core.atom     import Atom

class Lattice(object):
    def __init__(self, pars, sym):
        """
        Lattice object

        :param pars: [a,b,c,alpha,beta,gamma] (Angstrom & degrees)
        :param sym: spacegroup symmetry as 'F d -3 m' or 227
        """
        self.sym      = sym
        self.atom     = []
        self.Symmetry = None
        self.diffract = None
        self.genSymmetry()

        self.alpha = np.deg2rad(90)
        self.beta  = np.deg2rad(90)
        self.gamma = np.deg2rad(90)
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
        :return: 3x3 ATS Tensor

        """
        if self.Symmetry == None:
            self.genSymmetry()
        return self.Symmetry.genTensor(r,q)

    def xrmsTensor(self, M):
        """
        XRMS tensor in spherical approximation (Hannon et al.)

        :param M: magnetic structure factor (vector)
        :return:  3x3 XRMS Tensor
        """
        Fm = 1j * np.array([[0, M[2], -M[1]], [-M[2], 0, M[0]], [M[1], -M[0], 0]])
        return Fm

    def addAtom(self, element, r):
        rPos, _ = self.Symmetry.genPos(r)
        self.atom.append(Atom(element,r,rPos))

    def qMult(self,q):
        qTemp = np.dot(self.Symmetry.symOp, q)
        qEqv  = np.vstack({tuple(row) for row in qTemp})
        m = len(qEqv)

        return m, qEqv

    def genSF(self, Q):
        """
        Generate atomic structure factor for each atom

        :param Q: wavevector
        :return F: structure factor array
        """
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

    def recipLattice(self):
        V = 2.0 * self.a *self.b * self.b * np.sqrt(
                1 - np.cos(self.alpha)**2 - np.cos(self.beta)**2 - np.cos(self.gamma)**2
                + 2 * np.cos(self.alpha) * np.cos(self.beta) * np.cos(self.gamma))

        self.aStar = 2.0 / V * np.pi * self.b * self.c * np.sin(self.alpha)
        self.bStar = 2.0 / V * np.pi * self.a * self.c * np.sin(self.beta)
        self.cStar = 2.0 / V * np.pi * self.a * self.b * np.sin(self.gamma)

        self.alphaStar = 0.5* np.pi # np.arccos()
        self.betaStar  = 0.5* np.pi #
        self.gammaStar = 0.5* np.pi #

        return [self.aStar,self.bStar,self.cStar,self.alphaStar,self.betaStar,self.gammaStar]

#    def changeOfBasis(self):


    def angle(self,hkl1,hkl2):
        orthoRecip = np.array([[self.a,0,0],[0,self.b,0],[0,0,self.c]])

        p1 = np.dot(orthoRecip,hkl1)