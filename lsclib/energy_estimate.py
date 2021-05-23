# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 09:50:37 2020

@author: smithd24
"""

import pvlib
import pandas as pd
import numpy as np
import modified_haydavies_badescu as mh
from coordinate_transformations import sph2cart
from coordinate_transformations import cart2sph
import rotations as rm
import math


# data = 'data/78206_33.45_-112.06_tmy-2019.csv' # Phoenix data file
data = 'data/1252757_42.65_-73.74_tmy-2019.csv' # Albany data file
[metadata, TMYData] = pvlib.iotools.read_psm3(data)

def forecast_lsc(TMYData, metadata, DC_size, geometric_gain, opt_eff,
                 spect_mismatch, surface_tilt = 90, surface_azimuth = 180,
                 surface_twist = 0, albedo = .2):   

    nominal_dc = DC_size*1000
    ghi = TMYData['GHI']
    datetime = TMYData.index
    temp_air = TMYData['Temperature']
    wind_speed = TMYData['Wind Speed']
    latitude = metadata['Latitude']
    longitude = metadata['Longitude']
    tz = metadata['Time Zone']
    altitude = metadata['Elevation']
    name = metadata['Location ID']
    
    Location = pvlib.location.Location(latitude, longitude, tz, altitude, name)
    solar_position = Location.get_solarposition(datetime)
    solar_zenith = solar_position['zenith']
    solar_azimuth = solar_position['azimuth']
    alb = albedo
    
    # calculate plane of array irradiance
    irradiance = pvlib.irradiance.erbs(ghi, solar_zenith, datetime)
    dhi = irradiance['dhi']
    dni = irradiance['dni']
    dni_extra = pvlib.irradiance.get_extra_radiation(datetime)
    POA = pvlib.irradiance.get_total_irradiance(surface_tilt, surface_azimuth,
                                                solar_zenith, solar_azimuth,
                                                dni, ghi, dhi, dni_extra,
                                                airmass=None, albedo = alb,
                                                surface_type=None,
                                                model = 'haydavies')
    POA['poa_iso'],POA['poa_circ'] = mh.haydavies(surface_tilt, surface_azimuth,
                                               dhi, dni, dni_extra,
                                               solar_zenith, solar_azimuth,
                                               projection_ratio=None)
    
    # import mc model outputs
    
    # define values used for beam and circumsolar components of light
    lsc_stats_opteff = opt_eff[0]
    lsc_stats_m = spect_mismatch[0]
    
    # define values used for isotropic and ground-reflected components of light
    opt_eff_iso = opt_eff[1]
    opt_eff_grnd = opt_eff[2]
    spect_mismatch_iso = spect_mismatch[1]
    spect_mismatch_grnd = spect_mismatch[2]
    lsc_diff_mod = {'opt_eff': opt_eff_iso.loc[surface_tilt][0],
                    'spect_mismatch': spect_mismatch_iso.loc[surface_tilt][0]} 
    lsc_grnd_mod = {'opt_eff': opt_eff_grnd.loc[surface_tilt][0],
                    'spect_mismatch': spect_mismatch_grnd.loc[surface_tilt][0]} 
        
    # create set of interpolated data for theta and phi
    lsc_stats_opteff.index = lsc_stats_opteff.index.to_series().apply(lambda x: np.round(float(x),0))
    lsc_stats_opteff.columns = lsc_stats_opteff.columns.to_series().apply(lambda x: np.round(float(x),0))
    lsc_stats_interp = pd.DataFrame(columns = np.linspace(0,180,181), index = np.linspace(0,90,91)) 
    lsc_stats_interp = lsc_stats_interp.combine_first(lsc_stats_opteff)
    lsc_stats_interp = lsc_stats_interp.apply(pd.to_numeric)
    lsc_stats_interp = lsc_stats_interp.interpolate(method = 'values', axis = 0)
    lsc_stats_opteff = lsc_stats_interp.interpolate(method = 'values', axis = 1)
    
    lsc_stats_m.index = lsc_stats_m.index.to_series().apply(lambda x: np.round(float(x),0))
    lsc_stats_m.columns = lsc_stats_m.columns.to_series().apply(lambda x: np.round(float(x),0))
    lsc_stats_interp = pd.DataFrame(columns = np.linspace(0,180,181), index = np.linspace(0,90,91)) 
    lsc_stats_interp = lsc_stats_interp.combine_first(lsc_stats_m)
    lsc_stats_interp = lsc_stats_interp.apply(pd.to_numeric)
    lsc_stats_interp = lsc_stats_interp.interpolate(method = 'values', axis = 0)
    lsc_stats_m = lsc_stats_interp.interpolate(method = 'values', axis = 1)
    
    # establish matrix for spectral mismatch and opteff vs solar angle
    lsc_params = pd.concat([solar_zenith,solar_azimuth], axis=1)
    lsc_params['theta_i'] = np.nan
    lsc_params['phi_i'] = np.nan
    lsc_params['opt_eff'] = np.nan
    lsc_params['spect_mismatch'] = np.nan
    
    # determine hourly power estimates for LSC
    for i in range(len(POA.index)):
        # convert to cartesian
        solar_zenith_rad = math.radians(solar_zenith[i])
        solar_azimuth_rad = math.radians(solar_azimuth[i])
        vect_s = sph2cart(solar_zenith_rad, solar_azimuth_rad)
        # rotate by azimuth, tilt, twist
        surface_azimuth_rad = math.radians(surface_azimuth)
        surface_tilt_rad = math.radians(surface_tilt)
        surface_twist_rad = math.radians(surface_twist)
        vect = rm.z(-surface_azimuth_rad, vect_s)
        vect = rm.y(-surface_tilt_rad, vect)
        vect_i = rm.z(-surface_twist_rad, vect)
        # convert to spherical and map
        theta_i, phi_i = cart2sph(vect_i)
        lsc_params['theta_i'][i] = round(math.degrees(theta_i),0)
        lsc_params['phi_i'][i] = 180 - abs(round(math.degrees(phi_i),0))  
        if lsc_params['theta_i'][i] > 90:
            lsc_params['opt_eff'][i] = 0
            lsc_params['spect_mismatch'][i] = 0
        else:    
            lsc_params['opt_eff'][i] = lsc_stats_opteff.loc[lsc_params['theta_i'][i],
                                                              lsc_params['phi_i'][i]]
            lsc_params['spect_mismatch'][i] = lsc_stats_m.loc[lsc_params['theta_i'][i],
                                                              lsc_params['phi_i'][i]]
    
    # calculate global_eff without m for temperature calculations
    POA['poa_global_conc_cell'] = (lsc_params['opt_eff'] *
                              (POA['poa_direct'] + POA['poa_circ']) +
                              lsc_diff_mod['opt_eff'] *
                              POA['poa_iso'] +
                              lsc_grnd_mod['opt_eff'] *
                              POA['poa_ground_diffuse']) * geometric_gain
    
    temp_cell = pvlib.temperature.pvsyst_cell(POA['poa_global_conc_cell'],
                                              temp_air,
                                              wind_speed)
    # establish direct component used by solar cell
    POA['poa_direct_eff_cell'] = (lsc_params['opt_eff'] *
                                  lsc_params['spect_mismatch'] *
                                  POA['poa_direct']) * geometric_gain
    # establish circumsolar component used by solar cell
    POA['poa_circ_eff_cell'] = (lsc_params['opt_eff'] *
                                lsc_params['spect_mismatch'] *
                                POA['poa_circ']) * geometric_gain
    # establish isotropic component used by solar cell 
    POA['poa_iso_eff_cell'] = (lsc_diff_mod['opt_eff'] *
                              lsc_diff_mod['spect_mismatch'] *
                              POA['poa_iso']) * geometric_gain
    # establish ground-reflected component used by solar cell 
    POA['poa_ground_eff_cell'] = (lsc_grnd_mod['opt_eff'] *
                                  lsc_grnd_mod['spect_mismatch'] *
                                  POA['poa_ground_diffuse']) * geometric_gain
    # sum diffuse component used by solar cell 
    POA['poa_diffuse_eff_cell'] = (POA['poa_circ_eff_cell'] +
                                   POA['poa_iso_eff_cell'])
    # sum global component used by solar cell 
    POA['poa_global_eff_cell'] = (POA['poa_direct_eff_cell'] +
                                  POA['poa_diffuse_eff_cell'] +
                                  POA['poa_ground_eff_cell'])
    
    # Module parameters
    nom_power = 300
    L = 1.640  # m
    W = 0.992  # m
    alpha_sc = .005  # A/C
    gamma_ref = .910 #-0.00400  # 1/K - maybe %?
    mu_gamma =-.00038  # 1/K - maybe %?
    I_L_ref = 9.9  # A
    I_o_ref = 1.8E-11  # A
    R_sh_ref = 1700  # ohm
    R_sh_0 = 7000  # ohm
    R_s = 0.185  # ohm
    cells_in_series = 60
    R_sh_exp = 5.5
    EgRef = 1.121
    irrad_ref = 1000
    temp_ref = 25
    
    # calculate array nominal energy at STC
    num_modules = (nominal_dc/nom_power)
    modules_area = (num_modules*L*W)
    nom_power_per_area = nom_power/(L*W)
    nom_efficiency = nom_power_per_area/1000
    effective_irradiance = sum(POA['poa_global_eff_cell'])
    nom_energy = effective_irradiance*modules_area*nom_efficiency
    
    # calculate Pmpp
    five_params = []
    for i in range(len(temp_cell.index)):
        pvsyst_params = pvlib.pvsystem.calcparams_pvsyst(POA['poa_global_eff_cell'][i],
                                                        temp_cell[i], alpha_sc,
                                                        gamma_ref, mu_gamma,
                                                        I_L_ref, I_o_ref,
                                                        R_sh_ref, R_sh_0,
                                                        R_s, cells_in_series,
                                                        R_sh_exp, EgRef,
                                                        irrad_ref, temp_ref)
        sd_output = pvlib.pvsystem.singlediode(pvsyst_params[0],
                                                pvsyst_params[1],
                                                pvsyst_params[2],
                                                pvsyst_params[3],
                                                pvsyst_params[4])
        five_params.append(pvsyst_params)
    
    five_params = pd.DataFrame(five_params, columns = ['photocurrent',
                                                      'saturation_current',
                                                      'resistance_series',
                                                      'resistance_shunt',
                                                      'nNsVth'])
    
    sd_output = pvlib.pvsystem.singlediode(five_params.photocurrent,
                                            five_params.saturation_current,
                                            five_params.resistance_series,
                                            five_params.resistance_shunt,
                                            five_params.nNsVth)
    
    mpp_energy = sum(sd_output.p_mp)*num_modules
    estimates_hourly = POA
    estimates_hourly['dc_power'] = np.nan
    for i in range(0,len(POA)):
        estimates_hourly['dc_power'][i] = sd_output.p_mp[i]/(L*W)

    ref_yield = sum(POA['poa_global'])/irrad_ref
    spec_yield = mpp_energy/nominal_dc
    perf_ratio = spec_yield/ref_yield
    mpp_energy_m2 = mpp_energy/(L*W)
    
    print('Global horizontal irradiation:', round(sum(ghi)/1000), 'kWh/m^2')
    print('Effective irradiation on collectors:', round(effective_irradiance/1000), 'kWh/m^2')
    print('Array nominal energy (at STC effic.):', round(nom_energy/1000000,1), 'MWh')
    print('Array virtual energy at MPP:', round(mpp_energy/1000000,1), 'MWh')
    print('\nReference Yield:', round(ref_yield,1), 'kWh/kW')
    print('Final Yield:', round(spec_yield,1), 'kWh/kW')
    print('Performance Ratio:', round(perf_ratio,3))
    
    print('Global Effective POA:', 
          round(sum(POA['poa_global_eff_cell'])/1000, 1), 'kWh/m^2')
    print('Global Direct POA:', 
          round(sum(POA['poa_direct_eff_cell'])/1000, 1), 'kWh/m^2')
    print('Diffuse Effective POA:', 
          round(sum(POA['poa_diffuse_eff_cell'])/1000, 1), 'kWh/m^2')
    print('Ground-reflected Effective POA:', 
          round(sum(POA['poa_ground_eff_cell'])/1000, 1), 'kWh/m^2')
    print('Global Concentrated POA:',
          round(sum(POA['poa_global_conc_cell'])/1000, 1), 'kWh/m^2')
    
    return nominal_dc, mpp_energy_m2


def forecast_solar_panel(TMYData, metadata, DC_size, surface_tilt = 90,
                        surface_azimuth = 180, albedo = .2):

    nominal_dc = DC_size*1000
    ghi = TMYData['GHI']
    datetime = TMYData.index
    temp_air = TMYData['Temperature']
    wind_speed = TMYData['Wind Speed']
    latitude = metadata['Latitude']
    longitude = metadata['Longitude']
    tz = metadata['Time Zone']
    altitude = metadata['Elevation']
    name = metadata['Location ID']
    
    Location = pvlib.location.Location(latitude, longitude, tz, altitude, name)
    solar_position = Location.get_solarposition(datetime)
    solar_zenith = solar_position['zenith']
    solar_azimuth = solar_position['azimuth']
    aoi = pvlib.irradiance.aoi(surface_tilt, surface_azimuth,
                               solar_zenith, solar_azimuth)
    alb = albedo
    
    # calculate plane of array irradiance
    irradiance = pvlib.irradiance.erbs(ghi, solar_zenith, datetime)
    dhi = irradiance['dhi']
    dni = irradiance['dni']
    
    dni_extra = pvlib.irradiance.get_extra_radiation(datetime)
    POA = pvlib.irradiance.get_total_irradiance(surface_tilt, surface_azimuth,
                                                solar_zenith, solar_azimuth,
                                                dni, ghi, dhi, dni_extra,
                                                airmass=None, albedo = alb,
                                                surface_type=None,
                                                model = 'haydavies')
    POA['poa_iso'],POA['poa_circ'] = mh.haydavies(surface_tilt, surface_azimuth,
                                               dhi, dni, dni_extra,
                                               solar_zenith, solar_azimuth,
                                               projection_ratio=None)
    
    # calculate effective irradiance
    # iam_beam = pvlib.iam.physical(aoi, n=1.526, K=4.0, L=0.002)
    # iam_diff = pvlib.iam.marion_diffuse('physical', surface_tilt)
    iam_beam = pvlib.iam.ashrae(aoi)
    iam_diff = pvlib.iam.marion_diffuse('ashrae', surface_tilt)
    iam_sky = iam_diff['sky']
    iam_ground = iam_diff['ground']
    POA['poa_direct_eff'] = iam_beam*POA['poa_direct']
    POA['poa_circ_eff'] = iam_beam*POA['poa_circ']
    POA['poa_iso_eff'] = iam_sky*POA['poa_iso']
    POA['poa_ground_eff'] = iam_ground*POA['poa_ground_diffuse']
    POA['poa_diffuse_eff'] = (POA['poa_circ_eff'] +
                             POA['poa_iso_eff'])
    POA['poa_global_eff'] = (POA['poa_direct_eff'] +
                             POA['poa_diffuse_eff'] +
                             POA['poa_ground_eff'])
    
    # Module parameters
    nom_power = 300
    L = 1.640  # m
    W = 0.992  # m
    alpha_sc = .005  # A/C
    gamma_ref = .910 #-0.00400  # 1/K - maybe %?
    mu_gamma =-.00038  # 1/K - maybe %?
    I_L_ref = 9.9  # A
    I_o_ref = 1.8E-11  # A
    R_sh_ref = 1700  # ohm
    R_sh_0 = 7000  # ohm
    R_s = 0.185  # ohm
    cells_in_series = 60
    R_sh_exp = 5.5
    EgRef = 1.121
    irrad_ref = 1000
    temp_ref = 25
    
    # calculate array nominal energy at STC
    num_modules = nominal_dc/nom_power
    modules_area = num_modules*L*W
    nom_power_per_area = nom_power/(L*W)
    nom_efficiency = nom_power_per_area/1000
    effective_irradiance = sum(POA['poa_global_eff'])
    nom_energy = effective_irradiance*modules_area*nom_efficiency
    
    # calculate cell temperature with PVsyst approach
    temp_cell = pvlib.temperature.pvsyst_cell(POA['poa_global_eff'], temp_air,
                                              wind_speed)
    
    # calculate Pmpp
    five_params = []
    for i in range(len(temp_cell.index)):
        pvsyst_params = pvlib.pvsystem.calcparams_pvsyst(POA['poa_global_eff'][i],
                                                       temp_cell[i], alpha_sc,
                                                       gamma_ref, mu_gamma,
                                                       I_L_ref, I_o_ref,
                                                       R_sh_ref, R_sh_0,
                                                       R_s, cells_in_series,
                                                       R_sh_exp, EgRef,
                                                       irrad_ref, temp_ref)
        sd_output = pvlib.pvsystem.singlediode(pvsyst_params[0],
                                               pvsyst_params[1],
                                               pvsyst_params[2],
                                               pvsyst_params[3],
                                               pvsyst_params[4])
        five_params.append(pvsyst_params)
    
    five_params = pd.DataFrame(five_params, columns = ['photocurrent',
                                                      'saturation_current',
                                                      'resistance_series',
                                                      'resistance_shunt',
                                                      'nNsVth'])
    
    sd_output = pvlib.pvsystem.singlediode(five_params.photocurrent,
                                            five_params.saturation_current,
                                            five_params.resistance_series,
                                            five_params.resistance_shunt,
                                            five_params.nNsVth)
    
    mpp_energy = sum(sd_output.p_mp)*num_modules
    estimates_hourly = POA
    estimates_hourly['dc_power'] = np.nan
    for i in range(0,len(POA)):
        estimates_hourly['dc_power'][i] = sd_output.p_mp[i]/(L*W)
    
    ref_yield = sum(POA['poa_global'])/irrad_ref
    spec_yield = mpp_energy/nominal_dc
    perf_ratio = spec_yield/ref_yield
    mpp_energy_m2 = mpp_energy/(L*W)
    
    print('Global horizontal irradiation:', round(sum(ghi)/1000), 'kWh/m^2')
    print('Effective irradiation on collectors:', round(effective_irradiance/1000), 'kWh/m^2')
    print('Array nominal energy (at STC effic.):', round(nom_energy/1000000,2), 'MWh')
    print('Array virtual energy at MPP:', round(mpp_energy/1000000,2), 'MWh')
    print('\nReference Yield:', round(ref_yield,1), 'kWh/kW')
    print('Final Yield:', round(spec_yield,1), 'kWh/kW')
    print('Performance Ratio:', round(perf_ratio,3))
    
    print('Global Effective POA:', 
          round(sum(POA['poa_global_eff'])/1000, 1), 'kWh/m^2')
    print('Direct Effective POA:', 
          round(sum(POA['poa_direct_eff'])/1000, 1), 'kWh/m^2')
    print('Diffuse Effective POA:', 
          round(sum(POA['poa_diffuse_eff'])/1000, 1), 'kWh/m^2')
    print('Ground-reflected Effective POA:', 
          round(sum(POA['poa_ground_eff'])/1000, 1), 'kWh/m^2')

    return nominal_dc, mpp_energy_m2