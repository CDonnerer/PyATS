import numpy as np
import pylab as pl
from calcTh import calcTh
from calcAbs import calcAbs
from rotTensor import *
from sigpiBasis import sigpiBasis

def plotAzi(F,Q,azir,lat,E,n):
    psi = np.arange(0,360,1)
    Fsp = np.zeros(shape=(len(psi),2,2), dtype='complex128')
    Iss, Isp, Ips, Ipp = (np.zeros(len(psi)) for _ in range(4))
    T = rotTensor(Q, azir, lat)
    diffT = np.dot(T,np.dot(F,T.T))
    th = calcTh(Q,lat,E)
    if n != 0:
        absorb, _, _, _ = calcAbs(lat, n, Q, E, azir, psi)
    else:
        absorb=1

    for i in range(0,len(psi)):
        Fsp[i,:,:] = sigpiBasis(diffT,psi[i]-180,th)
        Iss[i] = abs(Fsp[i,0,0] * np.conjugate(Fsp[i,0,0]))
        Isp[i] = abs(Fsp[i,1,0] * np.conjugate(Fsp[i,1,0]))
        Ips[i] = abs(Fsp[i,0,1] * np.conjugate(Fsp[i,0,1]))
        Ipp[i] = abs(Fsp[i,1,1] * np.conjugate(Fsp[i,1,1]))

    pl.plot(psi,absorb*Isp, '-b', label='$\sigma-\pi$')
    pl.plot(psi,absorb*Iss, '-r', label='$\sigma-\sigma$')
    pl.plot(psi, absorb * (Iss+Isp), '-k', label='Total')
    pl.title('Q = ' + str(Q) +', $\psi_0$ = ' + str(azir))
    pl.xlim([0,360])
    pl.xlabel('$\psi$ (deg)')
    pl.ylabel('Intensity (arb. u.)')
    pl.show()