from genOp import genOp
from pointSym import pointSym
import numpy as np

def genTensor(sym,r,Q):
    symOp, symTr = genOp(sym)

    anisoT = np.ones(shape=(3,3)) - np.identity(3)
    isoT   = np.identity(3)
    Fk = np.zeros(shape=(3,3), dtype='complex128')

    for i in range(0,len(symOp)):
        Fk += np.dot(symOp[i],np.dot(anisoT, symOp[i].T)) \
              * np.exp(2j * np.pi * np.dot(Q, (symTr[i] + np.dot(symOp[i], r))))

    chop(Fk)
    norm = Fk.sum()/2

    return  Fk / norm

def chop(a):
    tol=1e-5
    a.real[abs(a.real) < tol] = 0
    a.imag[abs(a.imag) < tol] = 0
    return a