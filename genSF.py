# function to generate the structure factor
import numpy as np

def genSF(rPos,Q):
    F = 0

    for i in range(0, len(rPos)):
        F += np.exp(2j * np.pi * np.dot(Q, rPos[i].T))

    return F