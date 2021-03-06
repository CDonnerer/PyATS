import warnings
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
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
        Computes tensor that projects crystal frame onto diffraction frame

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

    def dspacing(self,Q):
        """
        Evaluate d-spacing of reflection

        :param Q: wavevector
        :return: d-spacing
        """

        if(len(np.ravel(Q)) == 3):
            Q = np.asarray(Q).reshape(1,3)

        #TODO generalise to triclinc symmetry
        gg = (Q[:,0] / self.lat.a)**2 + (Q[:,1] / self.lat.b)**2 + (Q[:,2] / self.lat.c)**2
        d = np.sqrt(1 / gg)
        return d

    def th(self, Q):
        """
        Calculate theta of reflection Q

        :param Q: in [h,k,l]
        :return: theta in radians
        """
        d = self.dspacing(Q)

        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            th = np.arcsin(self.lam / (2.0*d))
        return th

    def tth(self,Q):
        """
        Calculate two-theta of reflection Q

        :param Q: in [h,k,l]
        :return: two-theta in radians
        """
        return 2.0*self.th(Q)

    def xrms_tensor(self, M):
        """
        Calculate XRMS tensor in spherical approximation

        :param M: magnetic structure factor vector
        :return: 3 x 3 scattering tensor
        """
        Fm = 1j * np.array([[0, M[2], -M[1]], [-M[2], 0, M[0]], [M[1], -M[0], 0]])
        return Fm

    # ----- Absorption corrections ----- #

    def calc_absorption(self, Q, psi):
        th = self.th(Q)
        delta = self.calc_delta(Q, psi)
        alpha = th - delta
        beta  = th + delta
        abs   = (1 + np.sin(alpha) / np.sin(beta)) ** -1

        return abs, alpha, beta

    def calc_delta(self, Q, psi):
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

    def jones_matrix(self,F,th,psi):
        sig = np.array([np.sin(psi), np.cos(psi), 0])
        sigp = np.array([np.sin(psi), np.cos(psi), 0])
        pi = np.array([np.sin(th) * np.cos(psi), -np.sin(th) * np.sin(psi), -np.cos(th)])
        pip = np.array([-np.sin(th) * np.cos(psi), np.sin(th) * np.sin(psi), -np.cos(th)])

        m = np.array([[np.dot(sig, F.dot(sigp)), np.dot(sig, F.dot(pip))],
                      [np.dot(pi,  F.dot(sigp)), np.dot(pi,  F.dot(pip))]])
        return m

    def azimuthal_scan(self, F, Q, plot=True, absorb=False):
        psi = np.arange(0,360,1)
        T = self.orientate(Q)             # crystal orienation in diffraction frame
        diffT = np.dot(T, np.dot(F, T.T)) # project tensor into diffraction frame
        th = self.th(Q)

        abs_correction = np.ones(len(psi))
        if(absorb):
            abs_correction, _ , _ = self.calc_absorption(Q, psi) # check this

        M = self.jones_matrix(diffT, th, np.deg2rad(psi))

        Fs = np.einsum('i,ijm,j->m', self.pol, M, [1,0])
        Fp = np.einsum('i,ijm,j->m', self.pol, M, [0,1])

        if(plot):
            plt.figure(figsize=(16, 8))

            plt.plot(psi, np.abs(Fs) ** 2, '-b')
            plt.plot(psi, np.abs(Fp) ** 2, '-r')
            plt.plot(psi, np.abs(Fp) ** 2 + np.abs(Fs) ** 2, '-k')

            plt.title('Q = ' + str(Q) + ', $\psi_0$ = ' + str(self.azir))
            plt.xlabel('$\psi$ (deg)')
            plt.ylabel('Intensity (arb. u.)')
            plt.legend(['$\sigma \sigma^\prime$','$\sigma \pi$','sum'])
            plt.show()

        df = pd.DataFrame(data={'azi': psi,
                           'Fs': Fs, 'Fp': Fp},
                          columns=['azi','Fs','Fp'])
        return df

    # ----- Powder diffraction ----- #
    # still experimental

    def generate_reflections(self,n):
        TOL = 1e-8

        # generate all possible reflections up to (n,n,n)
        h = np.linspace(n,-n,2*n+1,dtype=np.int16)
        qs = np.array(np.meshgrid(h, h, h)).T.reshape(-1, 3)
        qs = qs[(qs != 0).any(axis=1)] # exclude (0,0,0)
        dspacings = self.dspacing(qs)

        # dspacing cutoff
        d_min = 0.5 * self.lam
        qs = qs[dspacings > d_min]
        ds = dspacings[dspacings > d_min]
        ds_inverse = 1 / (2.0*ds)

        # structure factor
        sfs = self.lat.generate_structure_factor(np.abs(qs)) # gives nans without abs? check

        # form factors
        for i, a in enumerate(self.lat.atom):
            fofa = a.formfactor(ds_inverse)
            sfs[:,i] *= fofa
        sf = np.abs(np.sum(sfs,axis=1))**2

        is_allowed = np.abs(sf) > TOL
        tths = 2*np.arcsin(self.lam * ds_inverse[is_allowed])
        qs = qs[is_allowed]
        sf = sf[is_allowed]

        # Lorentz polarisation factor
        LP = (1 + np.cos(tths) ** 2) / (8 * np.sin(0.5 * tths)**2 * np.cos(0.5 * tths))
        sf *= LP

        # Only unique reflections
        tth, idx, m = np.unique(tths, return_index=True, return_counts=True)
        hkl = qs[idx]
        int = m*sf[idx]
        int *= 100.0/max(int)

        df = pd.DataFrame(data={'tth': np.rad2deg(tth),
                           'h': hkl[:, 0], 'k': hkl[:, 1], 'l': hkl[:, 2],
                           'm': m, 'int': int},
                          columns=['tth','h','k','l','m','int'])
        return df

    def powder(self, n=20):
        df = self.generate_reflections(20)

        tth = np.asarray(df.tth)
        int = np.asarray(df.int)

        x=np.arange(0,90,.05)
        y=np.zeros_like(x)

        for t,s in zip(tth, int):
            y += self.lorz(x,s,t,.05)

        plt.figure(figsize=(16, 8))
        plt.plot(x,y)
        plt.xlabel('2 theta (deg)')
        plt.ylabel('Intensity (arb. u.)')
        plt.show()

    def lorz(self,x,a,c,w):
        L = a*w**2 / ( (x-c)**2 + w**2 )
        return L