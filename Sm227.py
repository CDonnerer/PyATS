from core.lattice import Lattice
import numpy as np
import pylab as pl

Sm227 = Lattice([10.3],'F d -3 m')

r = [0,0,0]
q = [14,0,0]

Fats = Sm227.genTensor(r,q)
Fm   = Sm227.xrmsTensor([1,0,0])

energy       =  11217
polarization = [1,0]
surface      = [1,0,0]
aziref       = [0,1,1]

Sm227.diffraction(energy, polarization, surface, aziref)
#Sm227.diffract.aziScan(Fats+Fm, q, absorb=False)

print(np.rad2deg(Sm227.diffract.tth([7,7,7])))
print(np.rad2deg(Sm227.diffract.tth([7.5,7.5,7.5])))

print(np.rad2deg(Sm227.diffract.tth([14,0,0])))
print(np.rad2deg(Sm227.diffract.tth([14.5,0.5,0.5])))


# psid = np.arange(0, 360, 1)
# th = Sm227.diffract.th(q)
# psi = np.deg2rad(psid)
#
# sig = np.array([np.sin(psi), np.cos(psi), 0])
# sigp = np.array([np.sin(psi), np.cos(psi), 0])
# pi  = np.array([np.sin(th) * np.cos(psi), -np.sin(th) * np.sin(psi), -np.cos(th)])
# pip = np.array([-np.sin(th) * np.cos(psi), np.sin(th) * np.sin(psi), -np.cos(th)])
#
# M = np.array([[np.dot(sig,Fats.dot(sigp)), np.dot(sig,Fats.dot(pip))],
#               [np.dot(pi,Fats.dot(sig)), np.dot(pi,Fats.dot(pip))]])
#
# print(M.shape)
#
# pin  = [1,0]
# pout = [0,1]
#
# #np.dot(pin,M.dot(pout))
#
# F = np.einsum('i,ijm,j->m',pin,M,pout)

#a = np.einsum('mij,jk,mlk->mil', self.symOp, anisoT, self.symOp)

#Fsp = np.zeros(shape=(len(psi), 2, 2), dtype='complex128')

Sm227.diffract.aziScan(Fats+0.5*Fm, q, absorb=True)
