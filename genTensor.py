from genOp import genOp
from pointSym import pointSym
import numpy as np

def genTensor(sym,r,Q):
    symOp, symTr = genOp(sym)

    q = np.zeros(shape=(6, 3, 3), dtype='complex128')
    q[0,:,:] = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    q[1,:,:] = np.array([[0, 0, 0], [0, 1, 0], [0, 0, 0]])
    q[2,:,:] = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]])
    q[3,:,:] = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
    q[4,:,:] = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    q[5,:,:] = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]])

    anisoT = np.ones(shape=(3,3)) - np.identity(3)
    isoT   = np.identity(3)
    #Fk = np.zeros(shape=(3,3), dtype='complex128')

    for j in range(0,6):
        Fk = np.zeros(shape=(3, 3), dtype='complex128')
        for i in range(0,len(symOp)):
            Fk += np.dot(symOp[i],np.dot(q[j,:,:], symOp[i].T)) \
                * np.exp(2j * np.pi * np.dot(Q, (symTr[i] + np.dot(symOp[i], r))))
        chop(Fk)
        print(Fk)

    print('---- Aniso T ----')
    Fk = np.zeros(shape=(3, 3), dtype='complex128')
    for i in range(0, len(symOp)):
        Fk += np.dot(symOp[i], np.dot(anisoT, symOp[i].T)) \
              * np.exp(2j * np.pi * np.dot(Q, (symTr[i] + np.dot(symOp[i], r))))
    chop(Fk)
    norm = Fk.sum() / 2

    # TODO catch if norm is zero
    return  Fk / norm

def chop(a):
    tol=1e-5
    a.real[abs(a.real) < tol] = 0
    a.imag[abs(a.imag) < tol] = 0
    return a