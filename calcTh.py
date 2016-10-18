import warnings
import numpy as np

def calcTh(Q, lat, E):

    lam = 1.23984193e4/E
    gg = (Q[0]/lat.a)**2 + (Q[1]/lat.b)**2 + (Q[2]/lat.c)**2
    d = np.sqrt(1/gg)

    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            th = np.arcsin(lam / (2 * d))
        except:
            print('Reflection Q = (',Q[0],Q[1],Q[2],') is not accessible at this energy!')
            return None

    return th