#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 09:50:35 2016

@author: lukas
"""

import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import axes3d
from scipy.spatial import ConvexHull, Delaunay
from matplotlib.path import Path
from scipy.stats import norm

def in_sphere(center, radius, x, y, z):

    square_dist = ((center[0] - x) ** 2 + (center[1] - y) ** 2 + (center[2] - z) ** 2 )
    return square_dist <= radius ** 2
    
def in_shell(center, radius_min, radius_max, point):
    x = point[0]
    y = point[1]
    z = point[2]
    square_dist = ((center[0] - x) ** 2 + (center[1] - y) ** 2 + (center[2] - z) ** 2 )
    return ((square_dist <= radius_max ** 2) and (square_dist >= radius_min ** 2))

def in_hull(p, hull):
    #http://stackoverflow.com/questions/16750618/whats-an-efficient-way-to-find-if-a-point-lies-in-the-convex-hull-of-a-point-cl
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull.points)
    return hull.find_simplex(p)>=0
        
def generate_rnd_spatial_points(n, minbound, maxbound):

    x_rnd = np.random.uniform(minbound[0],maxbound[0],n)
    y_rnd = np.random.uniform(minbound[1],maxbound[1],n)
    z_rnd = np.random.uniform(minbound[2],maxbound[2],n)

    rnd_points = np.zeros((n,3))
    rnd_points[:,0] = x_rnd
    rnd_points[:,1] = y_rnd
    rnd_points[:,2] = z_rnd

    return rnd_points
        
def generate_precipitation(pos, center=(0,0,0), radius=3):
    
    mc_matrix = 0
    mc_prec = 1

    pos = np.c_[pos, np.ones(pos[:,0].size) ]   
    
    for i,p in enumerate(pos):
        
        inside = in_sphere(center, radius, pos[i,0], pos[i,1], pos[i,2])
        if inside:
            pos[i,3] = mc_prec
        else:
            pos[i,3] = mc_matrix

    return pos

def sample_prec(pos, size_of_volume):

    # python adaption of the precipitation sampling from https://github.com/peterfelfer/AtomProbeTutorials    
    # for size of volume take the smallest absolute of the hull
    precRad = 3
    precConc = 90
    interfacialExcessL = 25
    interfacialExcessR = 25
    segWidth = 1
    rmax = np.sqrt(3* size_of_volume**2) *1.1
    mc_matrix = 0
    mc_prec = 1
    
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
    mcPrec[isSol_indices] = mc_prec
    mcPrec[noSol_indices] = mc_matrix
    
    posfile = np.c_[pos, np.ones(pos[:,0].size) ]    
    posfile[:,3] = mcPrec
    
    return posfile
    
fig = pl.figure()
ax = fig.add_subplot(111,projection='3d')

# text file is only the vertices of the blender exported obj
apt = np.genfromtxt('/home/lukas/sampled_atom_probe_data/scaled_mesh.txt')

points = np.zeros((apt[:,1].size,3))
points[:,0] = apt[:,1]
points[:,1] = apt[:,2]
points[:,2] = apt[:,3]

convex_hull = ConvexHull(points)

atomic_density = 20 #atoms/nm**3
# as the mesh is symmetrical take the abs of the smallest maxbounds
size_of_volume = np.min(convex_hull.max_bound)

number_of_atoms = int((size_of_volume * 2)**3 * atomic_density)

rnd_points = np.zeros((number_of_atoms,3))
rnd_points = generate_rnd_spatial_points(number_of_atoms, convex_hull.min_bound, convex_hull.max_bound)

posfile = np.zeros((number_of_atoms,3))

inside = np.where(in_hull(rnd_points, convex_hull))

posfile = rnd_points[inside]

mass_to_charge_matrix = 0
mass_to_charge_precipitation = 1

# create spherical precipitations
#posfile = sample_prec(posfile, size_of_volume)

posfile = generate_precipitation(posfile, center=(0,0,0), radius=8)

matrix = posfile[posfile[:,3] == 0]
prec = posfile[posfile[:,3] == 1]

n_matrix  = 10000
n_prec  = 100
ax.scatter(matrix[::n_matrix,0], matrix[::n_matrix,1], matrix[::n_matrix,2])
ax.scatter(prec[::n_prec,0], prec[::n_prec,1], prec[::n_prec,2], color='red')
for simplex in convex_hull.simplices:
    ax.plot(points[simplex, 0], points[simplex, 1],points[simplex, 2], 'k-')
pl.axis('equal')
ax.set_axis_off()
pl.show()

#np.savetxt('numpy_output_pos.txt',posfile)









