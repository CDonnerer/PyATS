from genTensor import *
from genPos import *
from rotTensor import *
from pointTensor import *
from lattice import Lattice
from plotAzi import plotAzi
from calcAbs import calcAbs
from qMult import *
from xrmsTensor import xrmsTensor
import numpy as np

sym  = 'P b c a'
ca214 = Lattice([5.38,5.630,11.722])
n = [0,0,1]

r    = [0,0,0]
q    = [0,1,3]
E    = 2967
azir = [0,1,0]
MS   = [0,1,0]

Ru_pos, _ = genPos(sym,r)

#print(Ru_pos)
#pointOp = pointSym(sym,r)
#print(pointOp)

Fats = genTensor(sym,r,q)
print(Fats)
Ftilt = np.array([[1,0,0],[0,0,1],[0,1,1]])


Fm = xrmsTensor(MS)

#FF = np.array([[0,0,0],[0,1,0],[0,0,0]])
#print(FF)

plotAzi(1.0*Fats+0.0*Fm+0.0*Ftilt,q,azir,ca214,E,0)