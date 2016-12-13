# ---- Eu2Ir2O7 powder diffraction ----- #

from core.lattice import Lattice
import numpy as np
import time
Eu227 = Lattice([10.3],'F d -3 m')

Eu227.addAtom('Eu',[0.5,0.5,0.5])
Eu227.addAtom('Ir',[0.0,0.0,0.0])
Eu227.addAtom('O' ,[0.375,0.375,0.125])
Eu227.addAtom('O' ,[0.332,0.125,0.125])

Eu227.diffraction(11217)

#np.rad2deg([1,2])

n = 5
qs  = []
ref = []
qs.append([0, 0, 0])
tol = 1e-5

t0 = time.time()
for h in range(1, n):
    for k in range(n):
        for l in range(n):
            try:
                tth = Eu227.diffract.tth([h, k, l])
                th = 0.5 * tth

                # Currently checking for equivalent reflections takes most time!
                # Is it needed? Could in theory loop over all h,k,l values?
                # For a list of powder reflections, equivalent ones need to be summed..

                # Best approach: Have list of symmetry-unique reflections available
                # before computing SFs.

                # m, qeq = Eu227.qMult([h, k, l])
                # truth = []
                # for r in qeq:
                # truth.append(any(np.sum(abs(np.asarray(r) - np.asarray(qs)), 1) < tol))
                if (True): #(not any(truth)):
                    qs.append([h, k, l])
                    #truth.clear()

                    SF = Eu227.genSF([h, k, l])
                    int = 0
                for i, f in enumerate(SF):
                    if (abs(f) < 1e-10):
                        continue
                    fofa = Eu227.atom[i].formfac(np.sin(th) / Eu227.diffract.lam)
                    int += f * fofa
                if (abs(int) > 1e-6):
                    LP = (1 + np.cos(tth) ** 2) / (8 * np.sin(th) ** 2 * np.cos(th))
                    ref.append([np.rad2deg(tth), [h, k, l], 1, 1 * LP * np.absolute(int) ** 2])
            except:
                continue
t1 = time.time()

print(1e3*(t1-t0))

Eu227.diffract.powder()

#ref = Eu227.diffract.genReflections(10)
#for r in ref:
#    print(r)
#Eu227.diffract.powder()


