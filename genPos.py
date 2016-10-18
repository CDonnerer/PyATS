import numpy as np
from genOp import genOp

def genPos(sym,r):
    if type(sym) is str:
        symOp, symTr = genOp(sym)
    else:
        symOp = sym[0]
        symTr = sym[1]

    rTemp = (np.dot(symOp,r) + symTr)%1
    rPos  = np.vstack({tuple(row) for row in rTemp})

    isMoved = np.sum(rTemp-r,axis=1) > 0

    return rPos, isMoved

#rPos = unique_rows(rTemp)
#def unique_rows(a):
#    a = np.ascontiguousarray(a)
#    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
#    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))