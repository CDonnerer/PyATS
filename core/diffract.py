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
        # Idea: Update local transformation matrix -> calculate everything from there

    # computes tensor that projects crystal frame onto diffraction frame
    def orientate(self, Q):
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
        psi = np.arctan2(-az[1], az[0]) # TODO this seems to off by 180 degrees. please check
        Ru3 = np.array([[np.cos(psi), -np.sin(psi), 0], [np.sin(psi), np.cos(psi), 0], [0, 0, 1]])

        return np.dot(Ru3, T)

    # generalise to triclinc symmetry
    def dspace(self,Q):
        gg = (Q[0] / self.lat.a)**2 + (Q[1] / self.lat.b)**2 + (Q[2] / self.lat.c)**2
        d = np.sqrt(1 / gg)
        return d

    def th(self, Q):
        d = self.dspace(Q)

        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            th = np.arcsin(self.lam / (2.0 * d))
        return th

    def tth(self,Q):
        return 2.0*self.th(Q)

    # XRMS tensor from magnetic SF
    def xrmsTensor(self, M):
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
        #npsi = a * psi[:, None]
        #n = diffVector(n, Q, azir, lat)  # n in lab frame at psi=0
        #   npsi = np.zeros(3)

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

    def aziScan(self, F, Q, absorb=False):
        psi = np.arange(0,360,1)
        Fsp = np.zeros(shape=(len(psi), 2, 2), dtype='complex128')
        Iss, Isp, = (np.zeros(len(psi)) for _ in range(2))
        T = self.orientate(Q)             # crystal orienation in diffraction frame
        diffT = np.dot(T, np.dot(F, T.T)) # project tensor into diffraction frame
        th = self.th(Q)

        absCorr = np.ones(len(psi))
        if(absorb):
            absCorr, _ , _ = self.calcAbs(Q, psi)

        for i, x in enumerate(psi):
            Fsp[i, :, :] = sigpiBasis(diffT, x - 180, th)
            #Iss[i] = absCorr[i]*abs(Fsp[i, 0, 0] * np.conjugate(Fsp[i, 0, 0]))
            Isp[i] = absCorr[i]*abs(Fsp[i, 1, 0] * np.conjugate(Fsp[i, 1, 0]))
            Iss[i] = 1*abs(Fsp[i, 1, 0] * np.conjugate(Fsp[i, 1, 0]))

        pl.plot(psi, Isp, '-b', label='$\sigma-\pi$')
        pl.plot(psi, Iss, '-r', label='$\sigma-\sigma$')
        pl.title('Q = ' + str(Q) +', $\psi_0$ = ' + str(self.azir))
        pl.xlim([0,360])
        pl.xlabel('$\psi$ (deg)')
        pl.ylabel('Intensity (arb. u.)')
        pl.show()

    # ----- Powder diffraction ----- #

    def genReflections(self,n):
        tol = 1e-8
        ref = []
        qs  = []
        qs.append([0,0,0])
        for h in range(1,n):
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
        reflex = self.genReflections(20)
        x=np.arange(0,90,.1)
        y=np.zeros_like(x)
        for r in reflex:
            y += self.lorz(x,r[3],r[0],.1)

        pl.plot(x,y)
        pl.xlabel('2 theta (deg)')
        pl.ylabel('Intensity (arb. u.)')
        pl.show()

    def lorz(self,x,a,c,w):
        L = a*w / ( (x-c)**2 + w**2 )
        return L