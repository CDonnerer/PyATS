import warnings
import pylab as pl
import numpy as np
from core.sigpiBasis import *
from core.lattice import *

class Diffract(object):
    def __init__(self, lattice, energy, polarization, surface, aziref, absorb):
        self.lat    = lattice
        self.e      = energy
        self.lam    = 1.23984193e4 / energy
        self.pol    = polarization
        self.n      = surface
        self.azir   = aziref
        self.absorb = absorb

    def orientate(self, Q):
        """
        computes tensor that projects crystal frame onto diffraction frame

        :param Q:
        :return:
        """
        Q = np.asarray(Q) / np.asarray([self.lat.a,self.lat.b,self.lat.c])
        Qnorm = Q / np.linalg.norm(Q)

        if (Q[0] == 0 and Q[2] == 0):
            zeta = 0
        else:
            zeta = np.arctan2(np.dot(Qnorm, np.array([0, 0, 1])), np.dot(Qnorm, np.array([1, 0, 0])))

        eta = np.arccos(np.dot(Qnorm, np.array([0, 1, 0])))
        T = np.array([[-np.cos(zeta) * np.cos(eta), np.sin(eta), -np.sin(zeta) * np.cos(eta)],
                      [ np.sin(zeta), 0, -np.cos(zeta)],
                      [-np.cos(zeta) * np.sin(eta),-np.cos(eta), -np.sin(zeta) * np.sin(eta)]])

        az = np.dot(T, self.azir)
        psi = np.arctan2(-az[1], az[0])
        Ru3 = np.array([[np.cos(psi), -np.sin(psi), 0], [np.sin(psi), np.cos(psi), 0], [0, 0, 1]])

        return np.dot(Ru3, T)

    def dspace(self,Q):
        """
        Evaluate d-spacing of reflection

        :param Q: wavevector
        :return: d-spacing
        """
        #TODO generalise to triclinc symmetry
        gg = (Q[0] / self.lat.a)**2 + (Q[1] / self.lat.b)**2 + (Q[2] / self.lat.c)**2
        d = np.sqrt(1 / gg)
        return d

    def th(self, Q):
        d = self.dspace(Q)

        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            th = np.arcsin(self.lam / (2.0*d))
        return th

    def tth(self,Q):
        return 2.0*self.th(Q)

    def xrmsTensor(self, M):
        """
        Calculate XRMS tensor in spherical approximation

        :param M: magnetic structure factor vector
        :return: 3 x 3 scattering tensor
        """
        Fm = 1j * np.array([[0, M[2], -M[1]], [-M[2], 0, M[0]], [M[1], -M[0], 0]])
        return Fm

    # ----- Absorption corrections ----- #

    def calcAbs(self, Q, psi):
        th = self.th(Q)
        delta = self.calcDelta(Q, psi)
        alpha = th - delta
        beta  = th + delta
        abs   = (1 + np.sin(alpha) / np.sin(beta)) ** -1

        return abs, alpha, beta

    def calcDelta(self, Q, psi):
        T = self.orientate(Q)   # crystal orienation in diffraction frame
        a = np.dot(T,self.n)
        delta = np.zeros(len(psi))

        for i in range(0, len(psi)):
            npsi = np.dot(self.rotZ(psi[i]), a)
            npro = self.nProj(npsi)
            delta[i] = np.arccos(np.dot(npro, np.array([0, 0, -1])))
            if npro[0] <= 0:
                delta[i] *= -1
        return delta

    def rotZ(self, angle):
        psi = np.deg2rad(angle)
        R = np.array([[np.cos(psi), -np.sin(psi), 0], [np.sin(psi), np.cos(psi), 0], [0, 0, 1]])
        return R

    def nProj(self, n):
        nplane = np.array([n[0], 0, n[2]])
        nPnorm = nplane / np.linalg.norm(nplane)
        return nPnorm

    # ----- Azi scan of tensor scattering ----- #

    def jonesMatrix(self,F,th,psi):
        sig = np.array([np.sin(psi), np.cos(psi), 0])
        sigp = np.array([np.sin(psi), np.cos(psi), 0])
        pi = np.array([np.sin(th) * np.cos(psi), -np.sin(th) * np.sin(psi), -np.cos(th)])
        pip = np.array([-np.sin(th) * np.cos(psi), np.sin(th) * np.sin(psi), -np.cos(th)])

        m = np.array([[np.dot(sig, F.dot(sigp)), np.dot(sig, F.dot(pip))],
                      [np.dot(pi,  F.dot(sigp)), np.dot(pi,  F.dot(pip))]])
        return m

    def aziScan(self, F, Q, absorb=False):
        psi = np.arange(0,360,1)
        T = self.orientate(Q)             # crystal orienation in diffraction frame
        diffT = np.dot(T, np.dot(F, T.T)) # project tensor into diffraction frame
        th = self.th(Q)

        absCorr = np.ones(len(psi))
        if(absorb):
            absCorr, _ , _ = self.calcAbs(Q, psi)

        M = self.jonesMatrix(diffT, th, np.deg2rad(psi))

        Fs = np.einsum('i,ijm,j->m', self.pol, M, [1,0])
        Fp = np.einsum('i,ijm,j->m', self.pol, M, [0,1])

        pl.plot(psi, np.abs(Fs)**2, '-b')
        pl.plot(psi, np.abs(Fp)**2, '-r')

        pl.title('Q = ' + str(Q) + ', $\psi_0$ = ' + str(self.azir))
        pl.xlabel('$\psi$ (deg)')
        pl.ylabel('Intensity (arb. u.)')
        pl.show()


    # ----- Powder diffraction ----- #
    # Experimental!

    def genReflections(self,n):
        tol = 1e-8
        ref = []
        qs  = []
        qs.append([0,0,0])
        for h in range(0,n):
            for k in range(n):
                for l in range(n):
                    try:
                        tth = self.tth([h,k,l])
                        th = 0.5 * tth
                        # m, qeq = self.lat.qMult([h, k, l])
                        # truth = []
                        # for r in qeq:
                        #    truth.append(any(np.sum(abs(np.asarray(r) - np.asarray(qs)), 1) < tol))
                        if(True): #not any(truth)):
                            qs.append([h,k,l])
                            #truth.clear()

                            SF = self.lat.genSF([h, k, l])
                            int = 0
                            for i, f in enumerate(SF):
                                if(abs(f) < 1e-10):
                                    continue
                                fofa = self.lat.atom[i].formfac(np.sin(th) / self.lam)
                                int += f * fofa
                            if (abs(int) > 1e-6):
                                LP = (1 + np.cos(tth) ** 2) / (8 * np.sin(th) ** 2 * np.cos(th))
                                ref.append([np.rad2deg(tth), [h, k, l], 1, 1 * LP * np.absolute(int) ** 2])
                    except:
                        continue
        return sorted(ref)

    def powder(self):
        reflex = self.genReflections(30)
        x=np.arange(0,120,.05)
        y=np.zeros_like(x)
        for r in reflex:
            y += self.lorz(x,r[3],r[0],.05)

        pl.plot(x,y)
        pl.xlabel('2 theta (deg)')
        pl.ylabel('Intensity (arb. u.)')
        pl.show()

    def lorz(self,x,a,c,w):
        L = a*w / ( (x-c)**2 + w**2 )
        return L