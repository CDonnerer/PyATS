from core.lattice import Lattice

Sm227 = Lattice([10.3],'F d -3 m')

r = [0,0,0]
q = [2,0,0]

Fats = Sm227.genTensor(r,q)
print(Fats)

energy       =  11217
polarization = [1,0]
surface      = [1,0,0]
aziref       = [0,0,1]


Sm227.diffraction(energy, polarization, surface, aziref)
#Sm227.diffract.aziScan(Fats,q)


# This is for later..
Sm227.addAtom([0,0,0],'Ir')
Sm227.diffract.powder()
