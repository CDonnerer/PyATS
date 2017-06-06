from core.lattice import Lattice
# ---- Pyrite ATS tensor ----- #

FeS2 = Lattice([5.417],'P a -3')

r = [0,0,0]
q = [0,1,1]
Fats = FeS2.genTensor(r,q)

print(Fats)