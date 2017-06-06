from core.lattice import Lattice
# ---- Eu2Ir2O7 powder diffraction ----- #

# create lattice
Eu227 = Lattice([10.3],'F d -3 m')

# atomic positions
Eu227.addAtom('Eu',[0.5,0.5,0.5])
Eu227.addAtom('Ir',[0.0,0.0,0.0])
Eu227.addAtom('O' ,[0.375,0.375,0.375])
Eu227.addAtom('O' ,[0.332,0.125,0.125])

# diffraction setup
Eu227.diffraction(11217)

# plot powder pattern
Eu227.diffract.powder()
