from rotTensor import *
import numpy as np

def calcAbs(lat,n,Q,E,azir,psi):
    th = calcTh(Q,lat,E)
    delta = calcDelta(lat, n, Q, azir, psi)
    alpha = th-delta
    beta  = th+delta
    abs = (1 + np.sin(alpha)/np.sin(beta))**-1

    return abs, alpha, beta, delta

def calcDelta( lat, n, Q, azir, psi):
    n = diffVector(n, Q, azir, lat)  # n in lab frame at psi=0
    npsi = np.zeros(3)
    delta = np.zeros(len(psi))

    for i in range(0, len(psi)):
        npsi = np.dot(rotZ(psi[i]), n)
        npro = nProj(npsi)
        delta[i] = np.arccos(np.dot(npro, np.array([0, 0, -1])))
        if npro[0] <= 0:
            delta[i] *= -1
    return delta

def rotZ(angle):
    psi = np.deg2rad(angle)
    R = np.array([[np.cos(psi), -np.sin(psi), 0], [np.sin(psi), np.cos(psi), 0], [0, 0, 1]])
    return R

def nProj(n):
    nplane = np.array([n[0], 0, n[2]])
    nPnorm = nplane / np.linalg.norm(nplane)
    return nPnorm