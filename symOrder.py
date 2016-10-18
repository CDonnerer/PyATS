# Determine the order of symmetry operator
import numpy as np

def symOrder(R, T):
    tol = 1e-5

    N = 1
    RN = R
    TN = T

    while (np.linalg.norm(RN - np.identity(3)) > tol or np.linalg.norm(TN) > tol) and N < 10:
        RN = np.dot(R, RN)
        TN = (R.dot(TN) + T) % 1
        N += 1

    return N