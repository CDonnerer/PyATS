from core.lattice import Lattice
# ---- Ca2RuO4 REXS ----- #

# create lattice and symmetry
Ca214 = Lattice([5.38,5.630,11.722],'P b c a')

r    = [0,0,0]
q    = [0,1,3]
MS   = [0,1,0]

# generate ATS and magnetic tensors
Fats = Ca214.genTensor(r,q)
Fm   = Ca214.xrmsTensor(MS)
print(Fats," ",Fm)

# setup diffraction experiment
n    = [0,0,1] # crystal surface
E    = 2967    # energy
pol  = [1,0]   # incident polarisation
azir = [0,1,0] # azimuthal reference

Ca214.diffraction(E,pol,n,azir)

# simulate azimuthal scan
Ca214.diffract.aziScan(Fats+0.5*Fm,q)