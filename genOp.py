# Calculate all symmetry operations of spacegroup
import numpy as np
from genSym import genSym
from symOrder import symOrder

def genOp(sym):
    tol = 1e-5
    gener, trans, symStr = genSym(sym)

    nGen = len(gener[1,1,:])
    symOp = []
    symTr = []
    symOp.append(np.identity(3))
    symTr.append(np.zeros(3))
    P = np.zeros(nGen, dtype=np.int)

    for i in range(0,nGen):
        R0 = gener[:,:,i]
        T0 = trans[:,i]
        P[i] = symOrder(R0,T0)

        R = np.identity(3)
        T = np.zeros(3)
        for j in range(0, P[i]-1):
            R = np.dot(R0,R)
            T = R0.dot(T) + T0
            nSym = len(symTr)

            for k in range(0,nSym):
                RS = np.dot(R,symOp[k])
                TS = (R.dot(symTr[k])+T)%1

                idxR = np.sum(np.sum(abs(symOp-RS),axis=1),axis=1) > tol
                idxT = np.sum(abs(symTr-TS),axis=1) > tol

                if all(idxT | idxR):
                    symTr.append(TS)
                    symOp.append(RS)
    return symOp, symTr