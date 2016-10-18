import numpy as np
from genOp import genOp
from genPos import genPos

def pointSym(sym,r):
    symOp, symTr = genOp(sym)

    _, isMoved = genPos(sym,r)

    opArray = np.asarray(symOp)
    pointOp = opArray[np.logical_not(isMoved)]

    return pointOp