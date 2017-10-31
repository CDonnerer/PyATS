from core.lattice import Lattice
import numpy as np
import pylab as pl

Sm227 = Lattice([10.3],'F d -3 m')

r = [0,0,0]
q = [2,0,0]

rTemp = (np.dot(Sm227.Symmetry.symOp, r) + Sm227.Symmetry.symTr) % 1

#T = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0.5]])


#b = np.exp(2j * np.pi * np.einsum('i,ji->j', q, (Sm227.Symmetry.symTr + np.dot(Sm227.Symmetry.symOp, r))))
#a = np.einsum('mij,jk,mlk->mil', Sm227.Symmetry.symOp, T, Sm227.Symmetry.symOp)
#Fk = np.einsum('ijk,i->jk', a, b)
#print(Fk)

b = np.exp(2j * np.pi * np.einsum('i,ji->j', q, (Sm227.Symmetry.symTr + np.dot(Sm227.Symmetry.symOp, r))))

print(np.exp(2j * np.pi * np.dot(q,[.0, 0, .0])))
print(np.exp(2j * np.pi * np.dot(q,[.0, 0.75, .75])))
print(np.exp(2j * np.pi * np.dot(q,[.75, 0, .75])))
print(np.exp(2j * np.pi * np.dot(q,[.75, .75, .0])))

#print(b[0:16:4])
#print(rTemp[0:16:3])



#(0, 0, 0) 2: (0, .75, .75)   3: (.75, 0, .75)   4: (.75, .75, 0)
#Fats = Sm227.genTensor(r,q)
#Fm   = Sm227.xrmsTensor([1,0,0])

#print(Fats)

#energy       =  11217
#polarization = [1,0]
#surface      = [1,0,0]
#aziref       = [0,1,1]

#Sm227.diffraction(energy, polarization, surface, aziref)

#print(np.rad2deg(Sm227.diffract.tth([7,7,7])))
#print(np.rad2deg(Sm227.diffract.tth([7.5,7.5,7.5])))

#print(np.rad2deg(Sm227.diffract.tth([14,0,0])))
#print(np.rad2deg(Sm227.diffract.tth([14.5,0.5,0.5])))

#Sm227.diffract.aziScan(Fats+0.5*Fm, q, absorb=True)
