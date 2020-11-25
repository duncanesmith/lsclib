# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 15:56:39 2019

@author: smithd24
"""

import math
import numpy as np
from lsclib import lsc_calcs as lc
from shapely.geometry import Polygon
from shapely.geometry import Point


def x(angle, vect): 
    """rotate a vector by an Euler angle with respect to x
    
    Parameters
    ----------
    angle : float
        Euler angle to rotate a vector with respect to x
    vect : list, NumPy array
        Coordinates encompassing the directional vector
    
    Returns
    -------
    vectx : NumPy array
        Rotated vector (after rotating vect by angle)
        
    Notes
    -----
    This is very slow - is there a faster way?
    """
    # Eulerian rotational matrix associated with the x-direction
    r_x = np.matrix([[1,         0,               0],
                    [0,         math.cos(angle), -math.sin(angle)],
                    [0,         math.sin(angle), math.cos(angle)]
                    ])

    vect = np.transpose(np.matrix(vect))
    vectx = r_x * vect  # determine rotated vector
    vectx = np.transpose(np.array(vectx))
    vectx = vectx[0]
    
    return vectx


def y(angle, vect):
    """rotate a vector by an Euler angle with respect to y
    
    Parameters
    ----------
    angle : float
        Euler angle to rotate a vector with respect to y
    vect : list, NumPy array
        Coordinates encompassing the directional vector
    
    Returns
    -------
    vecty : NumPy array
        Rotated vector (after rotating vect by angle)
        
    Notes
    -----
    This is very slow - is there a faster way?
    """
    # Eulerian rotational matrix associated with the y-direction    
    r_y = np.matrix([[math.cos(angle),    0,      math.sin(angle)],
                    [0,                  1,      0],
                    [-math.sin(angle),   0,      math.cos(angle)]
                    ])
    vect = np.transpose(np.matrix(vect))
    vecty = r_y * vect  # determine rotated vector
    vecty = np.transpose(np.array(vecty))
    vecty = vecty[0]
    return vecty


def z(angle, vect):
    """rotate a vector by an Euler angle with respect to z
    
    Parameters
    ----------
    angle : float
        Euler angle to rotate a vector with respect to z
    vect : list, NumPy array
        Coordinates encompassing the directional vector
    
    Returns
    -------
    vectz : NumPy array
        Rotated vector (after rotating vect by angle)
        
    Notes
    -----
    This is very slow - is there a faster way?
    """
    # Eulerian rotational matrix associated with the z-direction    
    r_z = np.matrix([[math.cos(angle),    -math.sin(angle),    0],
                    [math.sin(angle),    math.cos(angle),     0],
                    [0,                  0,                   1]
                    ])

    vect = np.transpose(np.matrix(vect))
    vectz = r_z * vect  # determine rotated vector
    vectz = np.transpose(np.array(vectz))
    vectz = vectz[0]
    return vectz


def rot(angle, vect, n):
    """rotate a vector by an Euler angle with respect to the surface normal of
    a boundary. This will rotate with respect to the x-axis or the y-axis.
    
    Parameters
    ----------
    angle : float
        Euler angle to rotate a vector with respect to n
    vect : list, NumPy array
        Coordinates encompassing the directional vector
    n : NumPy array
        Normal vector of a boundary
    
    Returns
    -------
    vect : NumPy array
        Rotated vector (after rotating vect by angle relative to n)
        
    Notes
    -----
    This is not an elegant way of doing this. This works for surfaces that are
    tilted with respect to the x-axis, and will work for surfaces with a normal
    that is parallel to the y-axis, but will not allow for anything else. For
    the code to be fully generalizable, this function will need to be expanded.
    """
    
    xvect = np.array([1,0,0])
    cross_test = np.cross(n, xvect)
    cross_test = list(cross_test)
    if cross_test == [0,0,0]:
        vect = y(angle, vect)  
    else:
        vect = x(angle, vect)

    return vect


def rot_point(angle, point, n):
    """rotate point into 2D plane in order to determine if it exists within a
    polygon. The Shapely library uses 2D geometry, so this is done in order to
    use it effectively for intersection calculations.
    
    Parameters
    ----------
    angle : float
        Euler angle to rotate a vector with respect to n
    point : NumPy array
        Coordinates encompassing a point (i.e. a potential intersection point)
    n : NumPy array
        Normal vector of a boundary
    
    Returns
    -------
    point_2d : Shapely Point object
        Shapely Point object in 2D coordinates
        
    Notes
    -----
    This is not an elegant way of doing this. This works for surfaces that are
    tilted with respect to the x-axis, and will work for surfaces with a normal
    that is parallel to the y-axis, but will not allow for anything else. For
    the code to be fully generalizable, this function will need to be expanded.
    """
    
    xvect = np.array([1,0,0])
    frontbacktest = lc.incidence_angle(n,xvect)
    # if this is a front or back surface of the LSC, rotate with respect to y
    if frontbacktest == 0 or frontbacktest == math.pi:
        point_2d = rot_point_y(angle, point)
    # otherwise, rotate with respect to x
    else:
        point_2d = rot_point_x(angle, point)
    
    return point_2d


def rot_poly(angle, polygon, n):
    """rotate polygon into 2D plane in order to determine if a point exists
    within it. The Shapely library uses 2D geometry, so this is done in order 
    to use it effectively for intersection calculations.
    
    Parameters
    ----------
    angle : float
        Euler angle to rotate a vector with respect to n
    polygon : NumPy array
        Coordinates encompassing a polygon (i.e. a boundary)
    n : NumPy array
        Normal vector of a boundary
    
    Returns
    -------
    poly_2d : Shapely Polygon object
        Shapely Polygon object in 2D coordinates
        
    Notes
    -----
    This is not an elegant way of doing this. This works for surfaces that are
    tilted with respect to the x-axis, and will work for surfaces with a normal
    that is parallel to the y-axis, but will not allow for anything else. For
    the code to be fully generalizable, this function will need to be expanded.
    """
       
    xvect = np.array([1,0,0])
    frontbacktest = lc.incidence_angle(n,xvect)
    # if this is a front or back surface of the LSC, rotate with respect to y    
    if frontbacktest == 0 or frontbacktest == math.pi:
        poly_2d = rot_poly_y(angle, polygon)
    # otherwise, rotate with respect to x    
    else:
        poly_2d = rot_poly_x(angle, polygon)
    
    return poly_2d


def rot_point_x(angle, point):
    """rotate point into 2D plane with respect to the x-direction.
    
    Parameters
    ----------
    angle : float
        Euler angle to rotate a vector with respect to n
    point : NumPy array
        Coordinates encompassing a point (i.e. a potential intersection point)

    Returns
    -------
    point_2d : Shapely Point object
        Shapely Point object in 2D coordinates
        
    Notes
    -----
    Deleting from an array is a relatively time consuming process I believe, it
    seems like this might be a place to target for increasing the speed of the
    code overall.
    """
       
    rotated_point = x(angle, point)  # rotate into xy plane
    point_2d = np.delete(rotated_point, 2, 0)  # delete z coordinate
    point_2d = Point(point_2d)  # convert point to shapely point

    return point_2d


def rot_poly_x(angle, polygon):
    """rotate polygon into 2D plane with respect to the x-direction.
    
    Parameters
    ----------
    angle : float
        Euler angle to rotate a vector with respect to n
    polygon : NumPy array
        Coordinates encompassing a polygon (i.e. a boundary)

    Returns
    -------
    poly_aligned : Shapely Polygon object
        Shapely Polygon object in 2D coordinates
        
    Notes
    -----
    Deleting from an array is a relatively time consuming process I believe, it
    seems like this might be a place to target for increasing the speed of the
    code overall.
    """
    polygon = np.array(polygon)
    poly_aligned = []
    # run rotations for all points in a polygon
    for i in range(0, len(polygon)):
        poly_aligned_3d = x(angle, polygon[i])  # rotate points into xy plane
        poly_aligned_2d = np.delete(poly_aligned_3d, 2, 0)  # delete z
        poly_aligned_2d = poly_aligned_2d.tolist()
        poly_aligned.append(poly_aligned_2d)
    poly_aligned = Polygon(poly_aligned)

    return poly_aligned


def rot_point_y(angle, point):
    """rotate point into 2D plane with respect to the y-direction.
    
    Parameters
    ----------
    angle : float
        Euler angle to rotate a vector with respect to n
    point : NumPy array
        Coordinates encompassing a point (i.e. a potential intersection point)

    Returns
    -------
    point_2d : Shapely Point object
        Shapely Point object in 2D coordinates
        
    Notes
    -----
    Deleting from an array is a relatively time consuming process I believe, it
    seems like this might be a place to target for increasing the speed of the
    code overall.
    """
    
    rotated_point = y(angle, point)  # rotate into xy plane
    point_2d = np.delete(rotated_point, 2, 0)  # delete z coordinate
    point_2d = Point(point_2d)  # convert point to shapely point

    return point_2d


def rot_poly_y(angle, polygon):
    """rotate polygon into 2D plane with respect to the y-direction.
    
    Parameters
    ----------
    angle : float
        Euler angle to rotate a vector with respect to n
    polygon : NumPy array
        Coordinates encompassing a polygon (i.e. a boundary)

    Returns
    -------
    poly_aligned : Shapely Polygon object
        Shapely Polygon object in 2D coordinates
        
    Notes
    -----
    Deleting from an array is a relatively time consuming process I believe, it
    seems like this might be a place to target for increasing the speed of the
    code overall.
    """
    polygon = np.array(polygon)
    poly_aligned = []
    # run rotations for all points in a polygon
    for i in range(0, len(polygon)):
        poly_aligned_3d = y(angle, polygon[i])  # rotate points into xy plane
        poly_aligned_2d = np.delete(poly_aligned_3d, 2, 0)  # delete z
        poly_aligned_2d = poly_aligned_2d.tolist()
        poly_aligned.append(poly_aligned_2d)
    poly_aligned = Polygon(poly_aligned)

    return poly_aligned