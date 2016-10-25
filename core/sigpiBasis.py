# Convert to sigma-pi 2x2 Jones matrix basis
# Assumes tensors is orientated in diffraction frame
# Returns amplitudes - square of intensities
# This requires wavelength and Q
# psi input is in degree

import numpy as np

def sigpiBasis(diffTens, psic, thBragg):
    psi = np.deg2rad(psic)

    pi  = np.array([ np.sin(thBragg)*np.cos(psi),-np.sin(thBragg)*np.sin(psi), -np.cos(thBragg)])
    pip = np.array([-np.sin(thBragg)*np.cos(psi), np.sin(thBragg)*np.sin(psi), -np.cos(thBragg)])

    sig  = np.array([np.sin(psi), np.cos(psi), 0])
    sigp = np.array([np.sin(psi), np.cos(psi), 0])

    sp = np.array([[ np.dot(sig,diffTens.dot(sigp)), np.dot(sig,diffTens.dot(pip)) ],
                   [ np.dot(pi ,diffTens.dot(sigp)), np.dot(pi ,diffTens.dot(pip)) ]])

    return sp