# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 12:26:43 2019

@author: smithd24
"""

import math


def sph2cart(theta, phi, r=1):
    """Converts spherical coords to cartesian"""

    x = r * math.sin(theta) * math.cos(phi)
    y = r * math.sin(theta) * math.sin(phi)
    z = r * math.cos(theta)

    vect = [x, y, z]

    return vect


def cart2sph(vect):
    """Converts cartesian coords to spherical"""
    
    theta = math.atan2(math.sqrt(vect[0]**2 + vect[1]**2), vect[2])
    phi = math.atan2(vect[1], vect[0])


    return [theta, phi]
