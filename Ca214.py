from core.lattice import Lattice


Ca214 = Lattice([5.38,5.630,11.722],'P b c a')

r    = [0,0,0]
q    = [0,1,3]
MS   = [0,1,0]

Fats = Ca214.genTensor(r,q)
Fm   = Ca214.xrmsTensor(MS)
print(Fats," ",Fm)

n    = [0,0,1]
E    =  2967
pol  = [1,0]
azir = [0,1,0]

Ca214.diffraction(E,pol,n,azir)
Ca214.diffract.aziScan(Fats+0.5*Fm,q)