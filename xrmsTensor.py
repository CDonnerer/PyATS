# Calculate the xrms tensor using magnetic structure factor,
# defined in the laboratory frame.
import numpy as np

def xrmsTensor(M):

    Fm = 1j*np.array([[0, M[2], -M[1]], [-M[2], 0, M[0]], [M[1],-M[0], 0]])
    return Fm

