from core.lattice import Lattice

Ca227 = Lattice([7.21,10.1,7.39,90,90,90], 'C m c m')

n = [1,1,1]

r    = [0,0,0]
q    = [0,0,5]
E    = 11217
azir = [0,1,0]
MS   = [0,0,1]


Fats = Ca227.genTensor(r, q)

#genTensor(sym,r,q)


#print(Fats)
#plotAzi(Fats,q,azir,sm227,E,0)

# * Make sure call to genOp is minimized, keep operators in mem.
# * Generalise rotTensor for hexagonal, monoclinic and triclinic -
#   this requires transforming to an orthogonal basis first
# * Restructure OO?