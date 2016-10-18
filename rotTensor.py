# Rotations tensor to transform from crystal to lab frame
# inputs are Q and azir (azimuthal reference)
# Set orientations such that azir lies in the scattering plane (outgoing)
# Polarisation vectors can then easily be defined.
# This currently only works for orthorhombic lattices!

import numpy as np

def rotTensor(Q,azir,lat):
    Q = np.asarray(Q) / np.array([lat.a,lat.b,lat.c])
    Qnorm = Q / np.linalg.norm(Q)

    if (Q[0]==0 and Q[2]==0):
        zeta = 0
    else:
        zeta = np.arctan2(np.dot(Qnorm, np.array([0, 0, 1])),np.dot(Qnorm, np.array([1, 0, 0])))

    eta = np.arccos( np.dot(Qnorm,np.array([0,1,0])) )
    T = np.array([ [-np.cos(zeta)*np.cos(eta), np.sin(eta),-np.sin(zeta)*np.cos(eta)],
                   [ np.sin(zeta),0, -np.cos(zeta)],
                   [-np.cos(zeta)*np.sin(eta),-np.cos(eta),-np.sin(zeta)*np.sin(eta)]])

    az = np.dot(T,azir)
    psi = np.arctan2(-az[1],az[0])
    Ru3 = np.array([[np.cos(psi),-np.sin(psi),0],[np.sin(psi),np.cos(psi),0],[0,0,1]])

    return chop(np.dot(Ru3,T))

def diffTensor(F,Q,azir,lat):
    T = rotTensor(Q,azir,lat)
    return np.dot(T, np.dot(F, T.T))

def diffVector(v,Q,azir,lat):
    T = rotTensor(Q,azir,lat)
    return np.dot(T, v)

def chop(a):
    tol=1e-5
    a.real[abs(a.real) < tol] = 0
    return a