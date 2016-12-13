from core.lattice import Lattice
import numpy as np
import time

Sm227 = Lattice([10.3],'F d -3 m')

r = [0,0,0]
q = [6,8,8]
Fats = Sm227.genTensor(r,q)

print(Fats)
#print(1e3*(t1-t0))
#n = 20
#qs = []
#ref = []
#qs.append([0, 0, 0])
#tol = 1e-5

#energy       = 11217
#polarization = [1,0]
#surface      = [1,1,1]
#aziref       = [1,0,0]

#Sm227.diffraction(energy, polarization, surface, aziref)
#Sm227.diffract.aziScan(Fats,q,absorb=True)