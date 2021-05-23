# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 12:11:46 2021

@author: smithd24
"""
import pvlib
import numpy as np


def haydavies(surface_tilt, surface_azimuth, dhi, dni, dni_extra,
              solar_zenith=None, solar_azimuth=None, projection_ratio=None):
    '''
    This is a slight modification from pvlib to output both the circumsolar
    and sky diffuse components of POA. The code below is almost entirely drawn
    from pvlib: https://pvlib-python.readthedocs.io/en/stable/generated/pvlib.irradiance.haydavies.html
    
    Determine diffuse irradiance from the sky on a tilted surface using
    Hay & Davies' 1980 model

    .. math::
        I_{d} = DHI ( A R_b + (1 - A) (\frac{1 + \cos\beta}{2}) )

    Hay and Davies' 1980 model determines the diffuse irradiance from
    the sky (ground reflected irradiance is not included in this
    algorithm) on a tilted surface using the surface tilt angle, surface
    azimuth angle, diffuse horizontal irradiance, direct normal
    irradiance, extraterrestrial irradiance, sun zenith angle, and sun
    azimuth angle.

    Parameters
    ----------
    surface_tilt : numeric
        Surface tilt angles in decimal degrees. The tilt angle is
        defined as degrees from horizontal (e.g. surface facing up = 0,
        surface facing horizon = 90)

    surface_azimuth : numeric
        Surface azimuth angles in decimal degrees. The azimuth
        convention is defined as degrees east of north (e.g. North=0,
        South=180, East=90, West=270).

    dhi : numeric
        Diffuse horizontal irradiance in W/m^2.

    dni : numeric
        Direct normal irradiance in W/m^2.

    dni_extra : numeric
        Extraterrestrial normal irradiance in W/m^2.

    solar_zenith : None or numeric, default None
        Solar apparent (refraction-corrected) zenith angles in decimal
        degrees. Must supply ``solar_zenith`` and ``solar_azimuth`` or
        supply ``projection_ratio``.

    solar_azimuth : None or numeric, default None
        Solar azimuth angles in decimal degrees. Must supply
        ``solar_zenith`` and ``solar_azimuth`` or supply
        ``projection_ratio``.

    projection_ratio : None or numeric, default None
        Ratio of angle of incidence projection to solar zenith angle
        projection. Must supply ``solar_zenith`` and ``solar_azimuth``
        or supply ``projection_ratio``.

    Returns
    --------
    sky_diffuse : numeric
        The sky diffuse component of the solar radiation.

    References
    -----------
    .. [1] Loutzenhiser P.G. et. al. "Empirical validation of models to
       compute solar irradiance on inclined surfaces for building energy
       simulation" 2007, Solar Energy vol. 81. pp. 254-267

    .. [2] Hay, J.E., Davies, J.A., 1980. Calculations of the solar
       radiation incident on an inclined surface. In: Hay, J.E., Won, T.K.
       (Eds.), Proc. of First Canadian Solar Radiation Data Workshop, 59.
       Ministry of Supply and Services, Canada.
    '''

    # if necessary, calculate ratio of titled and horizontal beam irradiance
    if projection_ratio is None:
        cos_tt = pvlib.irradiance.aoi_projection(surface_tilt, surface_azimuth,
                                solar_zenith, solar_azimuth)
        cos_tt = np.maximum(cos_tt, 0)  # GH 526
        cos_solar_zenith = pvlib.tools.cosd(solar_zenith)
        Rb = cos_tt / np.maximum(cos_solar_zenith, 0.01745)  # GH 432
    else:
        Rb = projection_ratio

    # Anisotropy Index
    AI = dni / dni_extra

    # these are the () and [] sub-terms of the second term of eqn 7
    term1 = 1 - AI
    term2 = 0.25 * (3 + pvlib.tools.cosd(2*surface_tilt))

    # sky_diffuse = dhi * (AI * Rb + term1 * term2)
    # sky_diffuse = np.maximum(sky_diffuse, 0)
    
    # this is appended onto the haydavies model
    isotropic = dhi * (term1 * term2)
    isotropic = np.maximum(isotropic, 0)
    circumsolar = dhi * (AI * Rb)
    circumsolar = np.maximum(circumsolar, 0)

    return isotropic, circumsolar