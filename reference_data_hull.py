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

def in_sphere(center, radius, point):
    x = point[0]
    y = point[1]
    z = point[2]
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
        
def generate_precipitation(center, radius, number_of_atoms):

    minbound_prec = center - np.array((radius,radius,radius))
    maxbound_prec = center + np.array((radius,radius,radius))
    
    prec_cube = generate_rnd_spatial_points(number_of_atoms, minbound_prec, maxbound_prec)
    
    # workaround don't get t right now
    spherical_prec = []
    for i,p in enumerate(prec_cube):
        
        inside = in_sphere(center, radius, p)
        if inside:
            spherical_prec.append(i)
            
    prec_sphere = prec_cube[spherical_prec]
    return prec_sphere
    
def generate_uniform_sphere_dist(radius, n):

    # http://stats.stackexchange.com/questions/85488/how-to-generate-random-points-in-the-volume-of-a-sphere-with-uniform-nearest-nei
    
    u = np.random.uniform(0,1,n)

    mu, sigma = 0, 1 # mean and standard deviation - variance shoul be on so sqrt(1) = 1
    x1 = np.random.normal(mu, sigma, n)
    x2 = np.random.normal(mu, sigma, n)
    x3 = np.random.normal(mu, sigma, n)
    
    x_rnd = (radius* (u**(1./3.))) /(np.sqrt(x1**2 + x2**2 + x3**2)) * x1
    y_rnd = (radius* (u**(1./3.))) /(np.sqrt(x1**2 + x2**2 + x3**2)) * x2
    z_rnd = (radius* (u**(1./3.))) /(np.sqrt(x1**2 + x2**2 + x3**2)) * x3

    rnd_points = np.zeros((n,3))
    rnd_points[:,0] = x_rnd
    rnd_points[:,1] = y_rnd
    rnd_points[:,2] = z_rnd

    return rnd_points
    
    
def unique_rows(data):
    uniq = np.unique(data.view(data.dtype.descr * data.shape[1]))
    return uniq.view(data.dtype).reshape(-1, data.shape[1])

def sample_prec(pos):


    # python adaption of the precipitation sampling from https://github.com/peterfelfer/AtomProbeTutorials    
    sz = 20
    precRad = 10
    precConc = 25
    interfacialExcessL = 25
    interfacialExcessR = 50
    segWidth = 1
    atomDens = 1
    rmax = np.sqrt(3* sz**2) *1.1
    mc1 = 0
    mc2 = 1
    
    scaler = 1000
#    # synthesis of position data
#    numAtom = (sz * 2)**3 * atomDens / scaler
#    pos = np.random.rand(numAtom,3)
#    #http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
#    pos = unique_rows(pos)
#    pos = np.array(pos) * sz * 2
#    pos = pos - sz
#    pos = pos[pos[:,2].argsort(),]
    
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
    
    return posfile
    
fig = pl.figure()
ax = fig.add_subplot(111,projection='3d')

# text file is only the vertices of the blender exported obj
apt = np.genfromtxt('/home/lukas/sampled_atom_probe_data/apt_hull.txt')

x = apt[:,1]
y = apt[:,2]
z = apt[:,3]

points = np.zeros((x.size,3))
points[:,0] = x
points[:,1] = y
points[:,2] = z

convex_hull = ConvexHull(points)

n = 10000
rnd_points = np.zeros((n,3))
rnd_points = generate_rnd_spatial_points(n, convex_hull.min_bound, convex_hull.max_bound)

posfile = np.zeros((n,3))

inside = np.where(in_hull(rnd_points, convex_hull))

posfile = rnd_points[inside]

mass_to_charge_matrix = 0
mass_to_charge_precipitation = 1

posfile = np.c_[posfile, np.ones(posfile[:,0].size) ]    
posfile[:,3] = mass_to_charge_matrix

# create spherical precipitations
sizes_of_precs = [5,10,20]
number_of_precs = 3

posfile = sample_prec(posfile)

#http://stackoverflow.com/questions/5408276/sampling-uniformly-distributed-random-points-inside-a-spherical-volume

#Generate a set of points uniformly distributed within a cube, then discard the ones whose distance from the center exceeds the radius of the desired sphere.

#center_prec1 = np.array((1,0,0))
#radius_prec1 = 0.75
#number_of_atoms_prec1 = 10000
#prec1_sphere = generate_precipitation(center_prec1,radius_prec1, number_of_atoms_prec1)
#
##prec1_sphere = generate_uniform_sphere_dist(radius_prec1, number_of_atoms_prec1)
##prec1_sphere += (center_prec1/2)
#
matrix = posfile[posfile[:,3] == 0]
prec = posfile[posfile[:,3] == 1]
ax.scatter(matrix[:,0], matrix[:,1], matrix[:,2])
ax.scatter(prec[:,0], prec[:,1], prec[:,2], color='red')
for simplex in convex_hull.simplices:
    ax.plot(points[simplex, 0], points[simplex, 1],points[simplex, 2], 'k-')
pl.axis('equal')
ax.set_axis_off()
pl.show()














