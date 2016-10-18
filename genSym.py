# Return symmetry generators of spacegroup
import numpy as np

def genSym(sym):
    f = open('dat_files/symmetry.dat')
    for line in f:
        if sym in line:
            symStr = line[19:]
            break
    f.close()

    vNew  = np.zeros(3)
    symOp = np.zeros(shape=(3, 3, 30))
    symTr = np.zeros(shape=(3, 30))
    nNew  = 0
    nOp   = 0
    nSign = 1

    j = 0
    while (j < len(symStr)):
        if symStr[j] == ',':
            symOp[nNew, :, nOp] = vNew
            vNew *= 0
            nSign = 1
            nNew  = nNew % 2 + 1
        elif symStr[j] == ';':
            symOp[nNew, :, nOp] = vNew
            vNew *= 0
            nSign = 1
            nNew  = 0
            nOp  += 1
        elif symStr[j] == 'x':
            vNew[0] = nSign
        elif symStr[j] == 'y':
            vNew[1] = nSign
        elif symStr[j] == 'z':
            vNew[2] = nSign
        elif symStr[j] == '-':
            nSign = -1
        elif symStr[j] == '+':
            nSign = 1
        elif symStr[j] == '1' or symStr[j] == '2' or symStr[j] == '3':
            symTr[nNew, nOp] = float(symStr[j]) / float(symStr[j + 2])
            j += 2
        j += 1

    symOp[nNew, :, nOp] = vNew
    symOp = symOp[:, :, 0:nOp + 1]
    symTr = symTr[:, 0:nOp + 1]
    return symOp, symTr, symStr