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

def in_sphere(center, radius, point):
    x = point[0]
    y = point[1]
    z = point[2]
    square_dist = ((center[0] - x) ** 2 + (center[1] - y) ** 2 + (center[2] - z) ** 2 )
    return square_dist <= radius ** 2

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

n = 2000
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

#http://stackoverflow.com/questions/5408276/sampling-uniformly-distributed-random-points-inside-a-spherical-volume

#Generate a set of points uniformly distributed within a cube, then discard the ones whose distance from the center exceeds the radius of the desired sphere.
#shareeditflag
center_prec1 = np.array((3,0,0))
radius_prec1 = 0.5
minbound_prec1 = center_prec1 - np.array((radius_prec1,radius_prec1,radius_prec1))
maxbound_prec1 = center_prec1 + np.array((radius_prec1,radius_prec1,radius_prec1))
number_of_atoms_prec1 = 300
prec1_cube = generate_rnd_spatial_points(number_of_atoms_prec1, minbound_prec1, maxbound_prec1)

# workaround don't get t right now
spherical_prec = []
for i,p in enumerate(prec1_cube):
    
    inside = in_sphere(center_prec1, radius_prec1, p)
    if inside:
        spherical_prec.append(i)
        
prec1_sphere = prec1_cube[spherical_prec]

#ax.scatter(posfile[:,0], posfile[:,1], posfile[:,2])
ax.scatter(prec1_sphere[:,0], prec1_sphere[:,1], prec1_sphere[:,2], color='red')
for simplex in convex_hull.simplices:
    ax.plot(points[simplex, 0], points[simplex, 1],points[simplex, 2], 'k-')

pl.show()














