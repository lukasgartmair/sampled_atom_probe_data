#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 15:44:35 2016

@author: lukas
"""

import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import axes3d
from scipy.stats import norm

# python adaption of the precipitation sampling from https://github.com/peterfelfer/AtomProbeTutorials

def unique_rows(data):
    uniq = np.unique(data.view(data.dtype.descr * data.shape[1]))
    return uniq.view(data.dtype).reshape(-1, data.shape[1])

sz = 20
precRad = 5
precConc = 25
baseConc = 1
interfacialExcessL = 25
interfacialExcessR = 50
segWidth = 1
atomDens = 1
rmax = np.sqrt(3* sz**2) *1.1
mc1 = 0
mc2 = 1

# synthesis of position data
numAtom = (sz * 2)**3 * atomDens
pos = np.random.rand(numAtom,3)
#http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
pos = unique_rows(pos)
pos = np.array(pos) * sz * 2
pos = pos - sz
pos = pos[pos[:,2].argsort(),]

## creating reference distribution
sample = np.random.rand(pos[:,0].size,1);
r = np.arange(0,rmax+1,sample[0])
concL = np.ones(r.size)
concL[r <= precRad] = precConc
concR = concL

##normpdf (X,mu,sigma)
excess =  norm.pdf(r, precRad, segWidth)
#
concL = concL + excess * interfacialExcessL
concR = concR + excess * interfacialExcessR
#
## sampling
atomRad = np.sqrt(np.sum(pos**2,1))
isL = pos[:,0] < 0
#
precPropVal = np.zeros(isL.size)
#
for i,L in enumerate(isL):
    if L == True:    
        precPropVal[i] = np.interp(atomRad[i],r,concL/100) # probability of an atom to be solute
    else:
        precPropVal[i] = np.interp(atomRad[i], r,concR/100)
#
isSol_indices = np.where(sample[:,0] < precPropVal)
noSol_indices = np.where(sample[:,0] >= precPropVal)
#
mcPrec = np.arange(precPropVal.size)
#
mcPrec[isSol_indices] = mc2
mcPrec[noSol_indices] = mc1

posfile = np.c_[pos, np.ones(pos[:,0].size) ]    
posfile[:,3] = mcPrec

matrix = posfile[posfile[:,3] == 0]
prec = posfile[posfile[:,3] == 1]

fig = pl.figure()
ax = fig.add_subplot(111,projection='3d')
ax.scatter(prec[:,0], prec[:,1], prec[:,2], color='blue')
pl.axis('equal')
pl.show()
