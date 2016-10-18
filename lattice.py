import numpy as np

class Lattice(object):
    def __init__(self, pars):
        self.alpha = np.deg2rad(90)
        self.beta  = np.deg2rad(90)
        self.gamma = np.deg2rad(90)

        if len(pars)==1:
            self.a=pars[0]
            self.b=self.a
            self.c=self.a

        elif len(pars)==3:
            self.a = pars[0]
            self.b = pars[1]
            self.c = pars[2]

        else:
            self.a = pars[0]
            self.b = pars[1]
            self.c = pars[2]
            self.alpha = np.deg2rad(pars[3])
            self.beta  = np.deg2rad(pars[4])
            self.gamma = np.deg2rad(pars[5])

        self.V = self.a * self.b * self.c \
                * np.sqrt( 1-np.cos(self.alpha)**2-np.cos(self.beta)**2-np.cos(self.gamma)**2\
                + 2*np.cos(self.alpha)*np.cos(self.beta)*np.cos(self.gamma) )