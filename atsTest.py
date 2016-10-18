from genTensor import *
from genPos import *
from genSym import *
from rotTensor import *
from pointTensor import *
from lattice import Lattice
from plotAzi import plotAzi
from calcAbs import calcAbs
from calcTh import *
from qMult import *
from xrmsTensor import xrmsTensor
import numpy as np

sym  = 'F d -3 m'
sm227 = Lattice([10.3])
n = [1,1,1]

r    = [0,0,0]
q    = [1,1,1]
E    = 1241
azir = [0,1,0]
MS   = [0,0,1]

sOp,sTr,sstr = genSym(sym)
print(sstr)
print(sOp[:,:,5])
print(sTr[:,5])

Fats = genTensor(sym,r,q)
#Fm   = xrmsTensor(MS)

#th = calcTh(q, sm227, E)
#print(2*np.rad2deg(th))

# print(Fats)
# print(Fm)
# T = rotTensor(Q,azir,sm227)

#plotAzi(0*Fats+1.0*Fm,q,azir,sm227,E,n)
# * Make sure call to genOp is minimized, keep operators in mem.
# * Generalise rotTensor for hexagonal, monoclinic and triclinic -
#   this requires transforming to an orthogonal basis first
# * Restructure OO?