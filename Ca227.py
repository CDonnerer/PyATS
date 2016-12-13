from core.lattice import Lattice

Ca227 = Lattice([7.21,10.1,7.39,90,90,90], 'C m c m')

n = [1,1,1]

r    = [0.25,0.25,0.25]
q    = [1,9,0]
E    = 10872
azir = [1,0,0]
MS   = [0,0,1]

Fats = Ca227.genTensor(r, q)
print(Fats)

Ca227.diffraction(E, [1,0], n, azir)
Ca227.diffract.aziScan(Fats, q)

# * Make sure call to genOp is minimized, keep operators in mem.
# * Generalise rotTensor for hexagonal, monoclinic and triclinic -
# * this requires transforming to an orthogonal basis first
# * Restructure OO?