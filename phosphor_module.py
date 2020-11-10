# -*- coding: utf-8 -*-
"""
Created on Wed May 22 10:23:21 2019

@author: smithd24
"""
import random
import numpy as np
import math
import time


def absorption_spectrum(spectrum, wave_len, wave_len_min, wave_len_max):
    """
    Determines wavelength of emitted bundle based upon the phosphor absorption
    spectrum (approximates the spectrum of sunlight) The normalized intensity
    per wavelength is called from excel and a random cumulative intensity is
    used to select a wavelength.

    Input:
    spectrum = normalized absorption spectrum

    Output:
    wave_len = incident wavelength of bundle
    """
#    start_time = time.perf_counter()

    probability = 0  # initialize probability of absorption

    if wave_len >= wave_len_min and wave_len <= wave_len_max:
        probability = spectrum.__call__(wave_len)
    
#    end_time = time.perf_counter() - start_time
#    print(end_time)
    return probability


def emission_spectrum(spectrum, spectrum_max):
    """
    Determines wavelength of emitted bundle based upon the phosphor emission
    spectrum (approximates the spectrum of sunlight) The normalized intensity
    per wavelength is called from excel and a random cumulative intensity is
    used to select a wavelength.

    Input:
    spectrum = normalized emission spectrum

    Output:
    wave_len = incident wavelength of bundle
    """
#    start_time = time.perf_counter()
    
    wave_len = 0  # initialize wave_len
    # if wavelength is specified directly then use it (spectrum is wavelength)
    if type(spectrum) is float:
        wave_len = spectrum

    # if spectrum is specified then solve for wavelength
    else:
        wave_len = spectrum.__call__(random.uniform(0, spectrum_max))
    
#    end_time = time.perf_counter() - start_time
#    print(end_time)
    
    return wave_len


def quantum_efficiency(wave_len, wave_len_em, qe):

    reset = 0
    absorptiontest = qe * (wave_len / wave_len_em)

    if random.uniform(0, 1) > absorptiontest:
        reset = 1

    return reset


def pathlength(extinction):
    
    # compute distance a bundle can travel before interacting with a phosphor
    # particle. "extinction" is the spectral absorption coefficient.
    pathlength = (-1/extinction)*math.log(random.uniform(0,1))

    return pathlength
