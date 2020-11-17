# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 14:58:33 2020

@author: smithd24
"""


def IoR_Sellmeier(wave_len):

    b1 = 1.0093    # Sellmeier coefficient for silicone
    c1 = 13185     # Sellmeier coefficient for silicone
    IoR = 1.0 + (b1*wave_len**2)/(wave_len**2-c1)
    IoR = IoR**.5  # Calculated index of refraction at a given wavelength

    return IoR
