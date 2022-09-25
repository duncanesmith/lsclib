# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 11:48:05 2020

@author: smithd24
"""

import pandas as pd
import numpy as np
import pvlib
import math
import time
import run

def LSC_stats_export(simulations, excel_name = 'LSC_stats_export.xlsx'):
    
    start_time = time.time()
    LSC = run.wedge(simulations, Einc = 1000, results = 'matrix')
    elapsed_time = time.time() - start_time
    print(elapsed_time)
    
    LSC_stats = pd.DataFrame(columns = ['theta','phi','optical_efficiency',
                             'spectral_mismatch','Isc_i','Isc_cell'])
    
    for j in range(len(LSC.columns)):
        for i in range(len(LSC.index)):
            new_row = {'theta': i*15, 'phi': j*15,
                       'optical_efficiency': LSC.iloc[i,j].optical_efficiency,
                       'spectral_mismatch': LSC.iloc[i,j].m,
                       'Isc_i': LSC.iloc[i,j].Isc_i,
                       'Isc_cell': LSC.iloc[i,j].Isc_cell}
            LSC_stats = LSC_stats.append(new_row, ignore_index = True)
    
    writer = pd.ExcelWriter(excel_name)
    LSC_stats.to_excel(writer,'Sheet1')
    writer.save()

    return LSC           


def module_from_pan_files(name):
    modules = pd.read_csv('data/PAN_files.csv')
    modules = modules.set_index('Model')
    module = modules.loc[name]
    return module


def dc_energy_estimate(DC_size, surface_tilt, weather_data, module_data,
                    surface_azimuth = 180, albedo = .2):
    
    #'data/1252757_42.65_-73.74_tmy-2019.csv'
    #'Mono 300 Wp 60 cells'
    
    # reformat weather data
    nominal_dc = DC_size*1000
    [metadata, TMYData] = pvlib.iotools.read_psm3(weather_data)
    ghi = TMYData['GHI']
    datetime = TMYData.index
    temp_air = TMYData['Temperature']
    wind_speed = TMYData['Wind Speed']
    latitude = metadata['Latitude']
    longitude = metadata['Longitude']
    tz = metadata['Time Zone']
    altitude = metadata['Elevation']
    name = metadata['Location ID']
    
    # reformat module data
    nom_power = module_data['Nom. Power']
    L = module_data['Length']/1000  # m
    W = module_data['Width']/1000  # m
    alpha_sc = module_data['muIsc']/100 # A/C
    gamma_ref = module_data['Gamma'] # 1/K
    mu_gamma = module_data['mu_gamma'] # 1/K
    I_L_ref = module_data['Isc']  # A
    I_o_ref = module_data['I_o_ref']  # A
    R_sh_ref = module_data['R shunt']  # ohm
    R_sh_0 = module_data['Rsh(G=0)']  # ohm
    R_s = module_data['R serie model']  # ohm
    cells_in_series = module_data['Nb cells series']
    R_sh_exp = module_data['Exponential parameter']
    EgRef = module_data['EgRef']
    irrad_ref = module_data['Gref'] # W/m^2 
    temp_ref = module_data['Tref'] # C
    
    # calculate solar angles
    Location = pvlib.location.Location(latitude, longitude, tz, altitude, name)
    solar_position = Location.get_solarposition(datetime)
    solar_zenith = solar_position['zenith']
    solar_azimuth = solar_position['azimuth']
    aoi = pvlib.irradiance.aoi(surface_tilt, surface_azimuth,
                               solar_zenith, solar_azimuth)
    
    # check for albedo data in TMY file
    if np.isnan(sum(TMYData['Surface Albedo'])):
        alb = albedo
    else:
        alb = TMYData['Surface Albedo']
    
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
    
    # calculate cell temperature
    temp_cell = pvlib.temperature.pvsyst_cell(POA['poa_global'], temp_air,
                                              wind_speed)
    
    # calculate effective irradiance
    #iam = pvlib.iam.physical(aoi, n=1.526, K=4.0, L=0.002)
    iam = pvlib.iam.ashrae(aoi)
    POA['poa_direct_eff'] = iam*POA['poa_direct']
    POA['poa_global_eff'] = POA['poa_diffuse'] + POA['poa_direct_eff']
    
    # calculate array nominal energy at STC
    num_modules = nominal_dc/nom_power
    modules_area = num_modules*L*W
    nom_power_per_area = nom_power/(L*W)
    nom_efficiency = nom_power_per_area/irrad_ref
    effective_irradiance = sum(POA['poa_global_eff'])
    nom_energy = effective_irradiance*modules_area*nom_efficiency
    
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
    dc_power = sd_output.p_mp
    mpp_energy = sum(dc_power)*num_modules
    
    print('Global horizontal irradiation:', round(sum(ghi)/1000), 'kWh/m^2')
    print('Effective irradiation on collectors:', round(effective_irradiance/1000), 'kWh/m^2')
    print('Array nominal energy (at STC effic.):', round(nom_energy/1000000), 'MWh')
    print('Array virtual energy at MPP:', round(mpp_energy/1000000), 'MWh')
    
    return POA, dc_power