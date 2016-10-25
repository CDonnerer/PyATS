import pylab as pl

from core.lattice import *
from core.sigpiBasis import sigpiBasis

class Diffract(object):
    def __init__(self, lattice, energy, polarization, surface, aziref):
        self.lat  = lattice
        self.e    = energy
        self.lam  = 1.23984193e4 / energy
        self.pol  = polarization
        self.n    = surface
        self.azir = aziref

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
                      [-np.cos(zeta) * np.sin(eta), -np.cos(eta), -np.sin(zeta) * np.sin(eta)]])

        az = np.dot(T, self.azir)
        psi = np.arctan2(-az[1], az[0]) # TODO this seems to off by 180 degrees. please check
        Ru3 = np.array([[np.cos(psi), -np.sin(psi), 0], [np.sin(psi), np.cos(psi), 0], [0, 0, 1]])

        return np.dot(Ru3, T)

    # generalise to triclinc symmetry
    def dspace(self,Q):
        gg = (Q[0] / self.lat.a) ** 2 + (Q[1] / self.lat.b) ** 2 + (Q[2] / self.lat.c) ** 2
        d = np.sqrt(1 / gg)
        return d

    #
    def theta(self,Q):
        d = self.dspace(Q)

        with np.warnings.catch_warnings():
            np.warnings.filterwarnings('error')
            try:
                th = np.arcsin(self.lam / (2.0 * d))
            except:
                print('Reflection Q = (', Q[0], Q[1], Q[2], ') is not accessible at this energy!')
                return None
        return th

    def tth(self,Q):
        return 2.0*self.theta(Q)


    def aziScan(self, F, Q):
        psi = np.arange(0,360,1)
        Fsp = np.zeros(shape=(len(psi), 2, 2), dtype='complex128')
        Iss, Isp, = (np.zeros(len(psi)) for _ in range(2))
        T = self.orientate(Q)             # crystal orienation in diffraction frame
        diffT = np.dot(T, np.dot(F, T.T)) # project tensor into diffraction frame
        th = self.theta(Q)

        for i in range(0, len(psi)):
            Fsp[i, :, :] = sigpiBasis(diffT, psi[i] - 180, th)
            Iss[i] = abs(Fsp[i, 0, 0] * np.conjugate(Fsp[i, 0, 0]))
            Isp[i] = abs(Fsp[i, 1, 0] * np.conjugate(Fsp[i, 1, 0]))

        pl.plot(psi, Isp, '-b', label='$\sigma-\pi$')
        pl.plot(psi, Iss, '-r', label='$\sigma-\sigma$')
        pl.title('Q = ' + str(Q) +', $\psi_0$ = ' + str(self.azir))
        pl.xlim([0,360])
        pl.xlabel('$\psi$ (deg)')
        pl.ylabel('Intensity (arb. u.)')
        pl.show()


    # ------ Blue Sky Thinking ------ #

    def powder(self):
        ref = []
        for h in range(5):
            for k in range(5):
                for l in range(5):
                    f = self.lat.genSF([h,k,l])
                    if (abs(f)>1e-10):
                        ref.append([h,k,l])
        print(ref)
        return ref