from genOp import genOp
from pointSym import pointSym
from genPos import genPos
import numpy as np

def symTensor(sym,r,q):
    tol=1e-5
    symOp, symTr = genOp(sym)

    qTemp = np.dot(symOp,q)
    qEqv  = np.vstack({tuple(row) for row in qTemp})

    print(qEqv)
    #isDifferent = np.sum(qTemp - q, axis=1) < 0

    #print(isDifferent)

#    Q = np.array([[1,0,0],[0,0,0],[0,0,0]])
#    Q = np.array([[0,1,0],[1,0,0],[0,0,0]])
#    o = 0

#    anisoT = np.ones(shape=(3,3)) - np.identity(3)
#    Fk = np.zeros(shape=(3,3), dtype='complex128')

#    for i in range(0,len(symOp)):
#        G = np.exp(2j * np.pi * np.dot(q, (symTr[i] + np.dot(symOp[i], r))))
#        Fk += np.dot(symOp[i],np.dot(Q, symOp[i].T)) * G
#        o += G

#    print(G)
#    print(Q)
#    print(Fk)

#    idxT = np.sum(Fk - Q, axis=1) < tol
#
#    print(idxT)
#    if all(idxT):
#        f = 1
#    else:
#        f=0
#    print(f)


#    chop(Fk)
#    norm = Fk.sum()/2

    return 0 #Fk # / norm

def chop(a):
    tol=1e-5
    a.real[abs(a.real) < tol] = 0
    a.imag[abs(a.imag) < tol] = 0
    return a