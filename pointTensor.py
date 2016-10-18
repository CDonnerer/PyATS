from pointSym import pointSym
import numpy as np
# TODO: This needs some thought...

def pointTensor(sym,r):
    tol = 1e-5
    symOp = pointSym(sym,r)

    q  = np.zeros(shape=(6,3,3), dtype='complex128')
    qT = np.zeros(shape=(6,3,3), dtype='complex128')

    q[0,:,:] = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    q[1,:,:] = np.array([[0, 0, 0], [0, 1, 0], [0, 0, 0]])
    q[2,:,:] = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]])
    q[3,:,:] = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
    q[4,:,:] = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    q[5,:,:] = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]])

    for i in range(0,len(q)):
        for j in range(0, len(symOp)):
            qT[i] += np.dot(symOp[j], np.dot(q[i], symOp[j].T))
        print(q[i])
        print(qT[i])

    invTens = np.zeros(shape=(3,3), dtype='complex128')
#    invTens.append(q[0])

    for i in range(0,len(q)):
        if (qT[0]-qT[i]).all() < tol:
 #           print(qT[0]-qT[i])
            invTens += q[i]

 #   print(invTens)

    return q[0]