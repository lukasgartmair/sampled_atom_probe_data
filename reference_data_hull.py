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

def in_hull(p, hull):
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
        
def generate_rnd_spatial_points(n,hull_minbound, hull_maxbound):

    x_rnd = np.random.uniform(hull_minbound[0],hull_maxbound[0],n)
    y_rnd = np.random.uniform(hull_minbound[1],hull_maxbound[1],n)
    z_rnd = np.random.uniform(hull_minbound[2],hull_maxbound[2],n)

    rnd_points = np.zeros((n,3))
    rnd_points[:,0] = x_rnd
    rnd_points[:,1] = y_rnd
    rnd_points[:,2] = z_rnd

    return rnd_points
        
fig = pl.figure()
ax = fig.add_subplot(111,projection='3d')

apt = np.genfromtxt('/home/lukas/master_thesis/codes/apt_hull.txt')

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

final_points = np.zeros((n,3))

inside = np.where(in_hull(rnd_points, convex_hull))

final_points = rnd_points[inside]

ax.scatter(final_points[:,0], final_points[:,1], final_points[:,2])
for simplex in convex_hull.simplices:
    ax.plot(points[simplex, 0], points[simplex, 1],points[simplex, 2], 'k-')

pl.show()

