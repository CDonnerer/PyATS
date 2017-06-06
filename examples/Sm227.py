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

print(np.rad2deg(Sm227.diffract.tth([7,7,7])))
print(np.rad2deg(Sm227.diffract.tth([7.5,7.5,7.5])))

print(np.rad2deg(Sm227.diffract.tth([14,0,0])))
print(np.rad2deg(Sm227.diffract.tth([14.5,0.5,0.5])))

Sm227.diffract.aziScan(Fats+0.5*Fm, q, absorb=True)
