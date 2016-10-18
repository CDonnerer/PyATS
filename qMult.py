# Find multiplicity of given reflection
import numpy as np
from genOp import genOp

def qMult(sym,r,q):
    symOp, symTr = genOp(sym)

    qTemp = np.dot(symOp,q)
    qEqv  = np.vstack({tuple(row) for row in qTemp})
    m = len(qEqv)

    return m, qEqv