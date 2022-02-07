# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:14:00 2019

@author: smithd24
"""

import math
import time
import xls_read_interpolate as xlsread
import lsc_classes as lsccls
import cProfile, pstats
import external_equations as eqn
import pandas as pd

import math
import time
import xls_read_interpolate as xlsread
import lsc_classes as lsccls
import cProfile, pstats
import external_equations as eqn
import pandas as pd

def wedge(trials, Einc = 1000, light_form = 'direct', results = 'single',
          theta_o = .000001, phi_o = .000001, tilt = 90,
          input_params = {'height': .007, 'film_thickness': .001,
                          'short_side_height': .0059, 'length': .03,
                          'mirror_gap': .007, 'width': .022,
                          'film_concentration': 27.73}):

    """Set up geometry to run trials for a wedge-shaped LSC. Specify
    vertices that make up the LSC. Each vertice will belong to a boundary and
    each boundary will belong to a volume. Using the LSC classes that have
    been defined, attribute characteristics to each boundary and volume.

    Parameters
    ----------
    trials : int
        Indicate number of bundles that will be used for the program
    Einc : float
        Indicates incident irradiance on the LSC. Default is 1000 W/m2.
    light_form : string
        Determines distribution of light entering the LSC.
        'direct'  - light enters at a fixed angle
        'diffuse' - diffuse light incident from the sky
        'ground'  - ground-reflected light incident from the ground
    results : string
        Determines if one result or if a matrix of incidence angle
        combinations is desired
    theta_o : float, optional
        Initial polar incidence angle. Default is zero. This is the polar angle
        relative to the LSC normal.
    phi_o : float, optional
        Initial azimuthal incidence angle. Default is zero. This is the azimuth
        angle relative to the LSC normal.
    tilt : float
        Tilt of LSC up from horizontal. This is input to compute the
        distribution of diffuse and ground-reflected light (when applicable)
    input_params : dict
        Dictionary of main LSC inputs.
    
    Returns
    -------
    LSC : object
        LSC object is returned after program has run. The LSC object will have
        the optical efficiency, short circuit current (for a cell within an
        LSC and a bare solar cell), and spectral mismatch factor available as
        attributes among a variety of other simulation results.
    
    Notes
    -----
    Having this written as one large script might not be the best. There are
    a large variety of inputs, so, instead of inputting them all here, these
    could be inputs to particular functions that make up "lsc_main". This would
    also demonstrate that a user is not limited to only the configuration
    detailed in the main script shown here.
    
    Phosphor particle should be added before iterating through incidence angle
    combinations.
    
    Starting volume/boundary process shoud be improved.
    
    errorcounts should be replaced by a for loop that runs automatically
    """
    
    # Initialize wedge-shaped geometry. Coordinates of individual vertices are
    # determined before they are assigned to individual boundaries.
    Height = input_params['height']
    Film_height = input_params['film_thickness']
    # distance up from 0 to the bottom point of the short side of the wedge
    Short_side_pos = input_params['short_side_height']
    Top_length = input_params['length']
    Mirror_gap_length = input_params['mirror_gap']
    W = input_params['width']
    precision = 16 
    hypotenuse = math.sqrt(Short_side_pos**2 + Top_length**2)
    angle = math.acos(Top_length/hypotenuse)
    
    L0 = 0
    H0_0 = 0
    H0_1 = Film_height
    H0_2 = Height
    
    L1 = Mirror_gap_length*math.cos(angle)
    H1_0 = Mirror_gap_length*math.sin(angle)
    H1_1 = H1_0 + Film_height
    
    L2 = Top_length
    H2_0 = Short_side_pos
    H2_1 = H2_0 + Film_height
    H2_2 = Height
    
    L1 = round(L1, precision)
    H1_0 = round(H1_0, precision)
    H1_1 = round(H1_1, precision)
    H2_1 = round(H2_1, precision)
    
    L = Top_length
    H = Height
    GG = L/H  # geometric gain (should eventually calculate from total PV area)
    mc_params = input_params
    mc_params['Geometric Gain'] = GG
    mc_params['Light Form'] = light_form
    mc_params['Trials'] = trials
    df_params = pd.DataFrame(mc_params, index = [0])
    
    # read in various excel data tables
    [abs_matrix, EQE_pv, IQE_pv, emi_source,
     abs_particle, emi_particle] = xlsread.excel_read()
    EQE_pv, EQE_pv_max = xlsread.spline(EQE_pv)
    IQE_pv, IQE_pv_max = xlsread.spline(IQE_pv)
    emi_source, emi_source_max = xlsread.spline(emi_source)
    abs_particle, abs_particle_max = xlsread.spline(abs_particle)
    emi_particle, emi_particle_max = xlsread.spline(emi_particle)
    
    # establish particle characteristics
    wave_len_min = 270  # minimum wavelength that can be absorbed by a particle
    wave_len_max = 500  # maximum wavelength that can be absorbed by a particle
    qe = 0.75           # quantum efficiency of a particle
    poa = .2            # probability of particle absorption (exp. value)
    extinction = 152.9*input_params['film_concentration'] # kg/m^3
    
    # establish matrix characteristics
    IoR = eqn.IoR_Sellmeier    # set index of refraction as constant or eqn
    abs_matrix, abs_matrix_max = xlsread.spline(abs_matrix)
    wave_len_min_matrix = 229   # minimum wavelength absorbed by matrix
    wave_len_max_matrix = 1100  # maximum wavelength absorbed by matrix
    
    
    # establish solar cell characteristics
    wave_len_min_pv = 350   # minimum wavelength absorbed by pv
    wave_len_max_pv = 1100  # maximum wavelength absorbed by pv
    
    # if running a combination of many theta_o and phi_o
    if light_form == 'direct' and results == 'matrix':
        
        data = {'':  [0.001, 15, 30, 45, 60, 75],
                 0:  [0, 0, 0, 0, 0, 0],
                15:  [0, 0, 0, 0, 0, 0],
                30:  [0, 0, 0, 0, 0, 0],
                45:  [0, 0, 0, 0, 0, 0],
                60:  [0, 0, 0, 0, 0, 0],
                75:  [0, 0, 0, 0, 0, 0],
                89.999:  [0, 0, 0, 0, 0, 0],
                105: [0, 0, 0, 0, 0, 0],
                120: [0, 0, 0, 0, 0, 0],
                135: [0, 0, 0, 0, 0, 0],
                150: [0, 0, 0, 0, 0, 0],
                165: [0, 0, 0, 0, 0, 0],
                180: [0, 0, 0, 0, 0, 0],} 
        
        # Convert the dictionary into DataFrame 
        df = pd.DataFrame(data)
        df.set_index('', inplace = True)
        df_oe = pd.DataFrame(data)
        df_oe.set_index('', inplace = True)
        df_m = pd.DataFrame(data)
        df_m.set_index('', inplace = True)
        theta_loop_count = len(df.index)
        phi_loop_count = len(df.columns)
    
    # if running a combination of many theta_o and phi_o
    if light_form != 'direct' and results == 'matrix':
        
        data = {'tilt': [ 5, 10, 20, 30, 40, 50, 60, 70, 80, 90], 
                'lsc_objects': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]}
        
        df = pd.DataFrame(data)
        df.set_index('tilt', inplace = True)
        df_oe = pd.DataFrame(data)
        df_oe.set_index('tilt', inplace = True)
        df_m = pd.DataFrame(data)
        df_m.set_index('tilt', inplace = True)
        
        theta_loop_count = len(df.index)
        phi_loop_count = 1
    
    # if expecting just one combination of inputs
    if results == 'single':
        
        theta_loop_count = 1
        phi_loop_count = 1
    
    for j in range(phi_loop_count):
        for i in range(theta_loop_count):
            start_time = time.time()
            lsc = lsccls.LSC()  # initialize LSC class
            
            # add phosphor particle
            particle = lsccls.Particle(poa, extinction, wave_len_min,
                                       wave_len_max, qe, abs_particle,
                                       emi_particle, emi_particle_max)
            
            # define dimensions/characteristics of Volume 0 - mirror gap
                        
            # input boundaries by setting boundary points of each
            bdy0a = [[0, L0, H0_0], [0, L1, H1_0], [W, L1, H1_0], [W, L0, H0_0]]
            bdy0b = [[0, L1, H1_0], [0, L1, H1_1], [W, L1, H1_1], [W, L1, H1_0]]
            bdy0c = [[0, L0, H0_1], [0, L1, H1_1], [W, L1, H1_1], [W, L0, H0_1]]
            bdy0d = [[0, L0, H0_0], [0, L0, H0_1], [W, L0, H0_1], [W, L0, H0_0]]
            bdy0e = [[0, L0, H0_0], [0, L1, H1_0], [0, L1, H1_1], [0, L0, H0_1]]
            bdy0f = [[W, L0, H0_0], [W, L1, H1_0], [W, L1, H1_1], [W, L0, H0_1]]
            bdys0 = [bdy0a, bdy0b, bdy0c, bdy0d, bdy0e, bdy0f]
             
            lsc.vol_list.append(lsccls.AbsorbingVolume(bdys0, 0, lsc, IoR,
                                                       abs_matrix,
                                                       wave_len_min_matrix,
                                                       wave_len_max_matrix))
            # add bottom surface
            lsc[0].bdy_list.append(lsccls.OpaqueBoundary(bdys0[0], lsc[0],
                                                         'specular', .05))
            # add right interface with film
            lsc[0].bdy_list.append(lsccls.TransparentBoundary(bdys0[1],
                                                              lsc[0]))
            # add interface with rest of matrix
            lsc[0].bdy_list.append(lsccls.TransparentBoundary(bdys0[2],
                                                              lsc[0]))
            # add left solar cell
            lsc[0].bdy_list.append(lsccls.PVBoundary(bdys0[3], lsc[0],
                                                     'diffuse', EQE_pv ,
                                                     0, wave_len_min_pv,
                                                     wave_len_max_pv))
            # add front mirror
            lsc[0].bdy_list.append(lsccls.OpaqueBoundary(bdys0[4], lsc[0],
                                                         'specular', .05))
            # add back mirror
            lsc[0].bdy_list.append(lsccls.OpaqueBoundary(bdys0[5], lsc[0],
                                                         'specular', .05))
            
            # define dimensions/characteristics of Volume 1 - phosphor film
            
            # input boundaries by setting boundary points of each
            bdy1a = [[0, L1, H1_0], [0, L2, H2_0], [W, L2, H2_0], [W, L1, H1_0]]
            bdy1b = [[0, L2, H2_0], [0, L2, H2_1], [W, L2, H2_1], [W, L2, H2_0]]
            bdy1c = [[0, L1, H1_1], [0, L2, H2_1], [W, L2, H2_1], [W, L1, H1_1]]
            bdy1d = [[0, L1, H1_0], [0, L1, H1_1], [W, L1, H1_1], [W, L1, H1_0]]
            bdy1e = [[0, L1, H1_0], [0, L2, H2_0], [0, L2, H2_1], [0, L1, H1_1]]
            bdy1f = [[W, L1, H1_0], [W, L2, H2_0], [W, L2, H2_1], [W, L1, H1_1]]
            bdys1 = [bdy1a, bdy1b, bdy1c, bdy1d, bdy1e, bdy1f]
            
            lsc.vol_list.append(lsccls.ParticleVolume(bdys1, 1, lsc, IoR,
                                                      abs_matrix, particle,
                                                      wave_len_min_matrix,
                                                      wave_len_max_matrix))
            # add bottom surface
            lsc[1].bdy_list.append(lsccls.OpaqueBoundary(bdys1[0], lsc[1],
                                                         'specular', .05))
            # add right mirror
            lsc[1].bdy_list.append(lsccls.OpaqueBoundary(bdys1[1], lsc[1],
                                                         'specular', .05))
            # add top surface
            lsc[1].bdy_list.append(lsccls.TransparentBoundary(bdys1[2], lsc[1]))
            # add left interface with mirror gap
            lsc[1].bdy_list.append(lsccls.TransparentBoundary(bdys1[3], lsc[1]))
            # add front mirror
            lsc[1].bdy_list.append(lsccls.OpaqueBoundary(bdys1[4], lsc[1],
                                                         'specular', .05))
            # add back mirror
            lsc[1].bdy_list.append(lsccls.OpaqueBoundary(bdys1[5], lsc[1],
                                                         'specular', .05))
            
            # define dimensions/characteristics of Volume 2 - rest of matrix
            
            # input boundaries by setting boundary points of each
            bdy2a = [[0, L0, H0_1], [0, L2, H2_1], [W, L2, H2_1], [W, L0, H0_1]]
            bdy2b = [[0, L2, H2_1], [0, L2, H2_2], [W, L2, H2_2], [W, L2, H2_1]]
            bdy2c = [[0, L0, H0_2], [0, L2, H2_2], [W, L2, H2_2], [W, L0, H0_2]]
            bdy2d = [[0, L0, H0_1], [0, L0, H0_2], [W, L0, H0_2], [W, L0, H0_1]]
            bdy2e = [[0, L0, H0_1], [0, L2, H2_1], [0, L2, H2_2], [0, L0, H0_2]]
            bdy2f = [[W, L0, H0_1], [W, L2, H2_1], [W, L2, H2_2], [W, L0, H0_2]]
            bdys2 = [bdy2a, bdy2b, bdy2c, bdy2d, bdy2e, bdy2f]
            
            # define volume
            lsc.vol_list.append(lsccls.AbsorbingVolume(bdys2, 2, lsc, IoR,
                                                       abs_matrix,
                                                       wave_len_min_matrix,
                                                       wave_len_max_matrix))
            # add interface with mirror gap and phosphor film
            lsc[2].bdy_list.append(lsccls.TransparentBoundary(bdys2[0],
                                                              lsc[2]))
            # add right mirror
            lsc[2].bdy_list.append(lsccls.OpaqueBoundary(bdys2[1], lsc[2],
                                                         'specular', .05))
            # add top surface
            lsc[2].bdy_list.append(lsccls.TransparentBoundary(bdys2[2],
                                                              lsc[2]))
            # add solar cell
            lsc[2].bdy_list.append(
                    lsccls.PVBoundary(bdys2[3], lsc[2], 'diffuse', EQE_pv , 0,
                                      wave_len_min_pv, wave_len_max_pv))
            # add front mirror
            lsc[2].bdy_list.append(lsccls.OpaqueBoundary(bdys2[4], lsc[2],
                                                         'specular', .05))
            # add back mirror
            lsc[2].bdy_list.append(lsccls.OpaqueBoundary(bdys2[5], lsc[2],
                                                         'specular', .05))
             
            # Prepare data inputs for LSC simulation
            if light_form == 'direct' and results == 'matrix':
                theta_o = df.index[i]
                phi_o = df.columns[j]
                
            if light_form != 'direct' and results == 'matrix':
                tilt = df.index[i]
                
            lsc.matching_pairs()    
            I = Einc*math.cos(math.radians(theta_o))*(L*W)
            theta_o = math.radians(theta_o + 180)  # adjust theta to head down
            phi_o = math.radians(phi_o + 90)       # adjust phi
            tilt = math.radians(tilt)              # adjust tilt
            
            # Run LSC trials, determining fate of every bundle
            starting_vol = len(lsc) - 1
            starting_bdy = 2
            
            lsc.main(trials, L, W, H, light_form, theta_o,
                     phi_o, tilt, starting_vol, starting_bdy, I,
                     emi_source, emi_source_max, particle)
            
            # Process data outputs from all LSC trials

            # determine if all bundles are accounted for
            bundle_count = particle.bundles_absorbed
            for vol_num in range(len(lsc)):
                try:
                    bundle_count = bundle_count + lsc[vol_num].bundles_absorbed
                except:
                    pass
                for bound_num in range(len(lsc[vol_num])):
                    try:
                        bundle_count = bundle_count + lsc[vol_num][bound_num].bundles_absorbed
                    except:
                        pass
                    try:
                        bundle_count = bundle_count + lsc[vol_num][bound_num].bundles_reflected
                    except:
                        pass
                    try:
                        bundle_count = bundle_count + lsc[vol_num][bound_num].bundles_refracted
                    except:
                        pass

            if bundle_count / trials != 1:
                print("\nENERGY IS NOT CONSERVED!!!!!")

            if light_form == 'direct' and results == 'maftrix':
                df.iloc[i,j] = lsc
                df_oe.iloc[i,j] = lsc.optical_efficiency
                df_m.iloc[i,j] = lsc.m
            
            if light_form != 'direct' and results == 'matrix':
                df.iloc[i,0] = lsc
                df_oe.iloc[i,0] = lsc.optical_efficiency
                df_m.iloc[i,0] = lsc.m
    
    if results == 'matrix':
        df_params.to_csv('mc_params.csv')
        lsc_outputs = {'parameters': mc_params,'opt_eff': df_oe,
                       'spect_mismatch':df_m, 'object':df,}
        if light_form =='direct':
            df_oe.to_csv('opt_eff_beam_circ.csv')
            df_m.to_csv('spect_mismatch_beam_circ.csv')
        if light_form =='diffuse':
            df_oe.to_csv('opt_eff_iso.csv')
            df_m.to_csv('spect_mismatch_iso.csv')
        if light_form =='ground':
            df_oe.to_csv('opt_eff_grnd.csv')
            df_m.to_csv('spect_mismatch_grnd.csv')    
        
    else:
        print(time.time() - start_time)
        print(light_form)
        print(lsc.optical_efficiency)
        print(lsc.m)
        print(tilt)
        lsc_outputs = {'parameters': mc_params,'opt_eff': lsc.optical_efficiency,
                       'spect_mismatch':lsc.m, 'object':lsc}

    return lsc_outputs


def wedge2(trials, Einc=1000, light_form='direct', results='single',
          theta_o=.000001, phi_o=.000001, tilt=90,
          input_params={'height': .007, 'film_thickness': .001,
                        'short_side_height': .0059, 'length': .03,
                        'mirror_gap': .007, 'width': .022,
                        'film_concentration': 27.73}):
    """Set up geometry to run trials for a wedge-shaped LSC. Specify
    vertices that make up the LSC. Each vertice will belong to a boundary and
    each boundary will belong to a volume. Using the LSC classes that have
    been defined, attribute characteristics to each boundary and volume.

    Parameters
    ----------
    trials : int
        Indicate number of bundles that will be used for the program
    Einc : float
        Indicates incident irradiance on the LSC. Default is 1000 W/m2.
    light_form : string
        Determines distribution of light entering the LSC.
        'direct'  - light enters at a fixed angle
        'diffuse' - diffuse light incident from the sky
        'ground'  - ground-reflected light incident from the ground
    results : string
        Determines if one result or if a matrix of incidence angle
        combinations is desired
    theta_o : float, optional
        Initial polar incidence angle. Default is zero. This is the polar angle
        relative to the LSC normal.
    phi_o : float, optional
        Initial azimuthal incidence angle. Default is zero. This is the azimuth
        angle relative to the LSC normal.
    tilt : float
        Tilt of LSC up from horizontal. This is input to compute the
        distribution of diffuse and ground-reflected light (when applicable)
    input_params : dict
        Dictionary of main LSC inputs.

    Returns
    -------
    LSC : object
        LSC object is returned after program has run. The LSC object will have
        the optical efficiency, short circuit current (for a cell within an
        LSC and a bare solar cell), and spectral mismatch factor available as
        attributes among a variety of other simulation results.

    Notes
    -----
    Having this written as one large script might not be the best. There are
    a large variety of inputs, so, instead of inputting them all here, these
    could be inputs to particular functions that make up "lsc_main". This would
    also demonstrate that a user is not limited to only the configuration
    detailed in the main script shown here.

    Phosphor particle should be added before iterating through incidence angle
    combinations.

    Starting volume/boundary process shoud be improved.

    errorcounts should be replaced by a for loop that runs automatically
    """

    # Initialize wedge-shaped geometry. Coordinates of individual vertices are
    # determined before they are assigned to individual boundaries.
    Height = input_params['height']
    Film_height = input_params['film_thickness']
    # distance up from 0 to the bottom point of the short side of the wedge
    Short_side_pos = input_params['short_side_height']
    Top_length = input_params['length']
    Mirror_gap_length = input_params['mirror_gap']
    W = input_params['width']
    precision = 16
    hypotenuse = math.sqrt(Short_side_pos ** 2 + Top_length ** 2)
    angle = math.acos(Top_length / hypotenuse)

    L0 = 0
    H0_0 = 0
    H0_1 = Film_height
    H0_2 = Height

    L1 = Mirror_gap_length * math.cos(angle)
    H1_0 = Mirror_gap_length * math.sin(angle)
    H1_1 = H1_0 + Film_height

    L2 = Top_length
    H2_0 = Short_side_pos
    H2_1 = H2_0 + Film_height
    H2_2 = Height

    L1 = round(L1, precision)
    H1_0 = round(H1_0, precision)
    H1_1 = round(H1_1, precision)
    H2_1 = round(H2_1, precision)

    L = Top_length
    H = Height
    GG = L / H  # geometric gain (should eventually calculate from total PV area)
    mc_params = input_params
    mc_params['Geometric Gain'] = GG
    mc_params['Light Form'] = light_form
    mc_params['Trials'] = trials
    df_params = pd.DataFrame(mc_params, index=[0])

    # read in various excel data tables
    [abs_matrix, EQE_pv, IQE_pv, emi_source,
     abs_particle, emi_particle] = xlsread.excel_read()
    EQE_pv, EQE_pv_max = xlsread.spline(EQE_pv)
    IQE_pv, IQE_pv_max = xlsread.spline(IQE_pv)
    emi_source, emi_source_max = xlsread.spline(emi_source)
    abs_particle, abs_particle_max = xlsread.spline(abs_particle)
    emi_particle, emi_particle_max = xlsread.spline(emi_particle)

    # establish particle characteristics
    wave_len_min = 270  # minimum wavelength that can be absorbed by a particle
    wave_len_max = 500  # maximum wavelength that can be absorbed by a particle
    qe = 0.75  # quantum efficiency of a particle
    poa = .2  # probability of particle absorption (exp. value)
    extinction = 152.9 * input_params['film_concentration']  # kg/m^3

    # establish matrix characteristics
    IoR = eqn.IoR_Sellmeier  # set index of refraction as constant or eqn
    abs_matrix, abs_matrix_max = xlsread.spline(abs_matrix)
    wave_len_min_matrix = 229  # minimum wavelength absorbed by matrix
    wave_len_max_matrix = 1100  # maximum wavelength absorbed by matrix

    # establish solar cell characteristics
    wave_len_min_pv = 350  # minimum wavelength absorbed by pv
    wave_len_max_pv = 1100  # maximum wavelength absorbed by pv

    # if running a combination of many theta_o and phi_o
    if light_form == 'direct' and results == 'matrix':
        data = {'': [0.001, 15, 30, 45, 60, 75],
                0: [0, 0, 0, 0, 0, 0],
                15: [0, 0, 0, 0, 0, 0],
                30: [0, 0, 0, 0, 0, 0],
                45: [0, 0, 0, 0, 0, 0],
                60: [0, 0, 0, 0, 0, 0],
                75: [0, 0, 0, 0, 0, 0],
                89.999: [0, 0, 0, 0, 0, 0],
                105: [0, 0, 0, 0, 0, 0],
                120: [0, 0, 0, 0, 0, 0],
                135: [0, 0, 0, 0, 0, 0],
                150: [0, 0, 0, 0, 0, 0],
                165: [0, 0, 0, 0, 0, 0],
                180: [0, 0, 0, 0, 0, 0], }

        # Convert the dictionary into DataFrame
        df = pd.DataFrame(data)
        df.set_index('', inplace=True)
        df_oe = pd.DataFrame(data)
        df_oe.set_index('', inplace=True)
        df_m = pd.DataFrame(data)
        df_m.set_index('', inplace=True)
        theta_loop_count = len(df.index)
        phi_loop_count = len(df.columns)

    # if running a combination of many theta_o and phi_o
    if light_form != 'direct' and results == 'matrix':
        data = {'tilt': [5, 10, 20, 30, 40, 50, 60, 70, 80, 90],
                'lsc_objects': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]}

        df = pd.DataFrame(data)
        df.set_index('tilt', inplace=True)
        df_oe = pd.DataFrame(data)
        df_oe.set_index('tilt', inplace=True)
        df_m = pd.DataFrame(data)
        df_m.set_index('tilt', inplace=True)

        theta_loop_count = len(df.index)
        phi_loop_count = 1

    # if expecting just one combination of inputs
    if results == 'single':
        theta_loop_count = 1
        phi_loop_count = 1

    for j in range(phi_loop_count):
        for i in range(theta_loop_count):
            start_time = time.time()
            lsc = lsccls.LSC()  # initialize LSC class

            # add phosphor particle
            particle = lsccls.Particle(poa, extinction, wave_len_min,
                                       wave_len_max, qe, abs_particle,
                                       emi_particle, emi_particle_max)

            # define dimensions/characteristics of Volume 0 - mirror gap

            # input boundaries by setting boundary points of each
            bdy0a = [[0, L0, H0_0], [0, L1, H1_0], [W, L1, H1_0], [W, L0, H0_0]]
            bdy0b = [[0, L1, H1_0], [0, L1, H1_1], [W, L1, H1_1], [W, L1, H1_0]]
            bdy0c = [[0, L0, H0_1], [0, L1, H1_1], [W, L1, H1_1], [W, L0, H0_1]]
            bdy0d = [[0, L0, H0_0], [0, L0, H0_1], [W, L0, H0_1], [W, L0, H0_0]]
            bdy0e = [[0, L0, H0_0], [0, L1, H1_0], [0, L1, H1_1], [0, L0, H0_1]]
            bdy0f = [[W, L0, H0_0], [W, L1, H1_0], [W, L1, H1_1], [W, L0, H0_1]]
            bdys0 = [bdy0a, bdy0b, bdy0c, bdy0d, bdy0e, bdy0f]

            lsc.vol_list.append(lsccls.AbsorbingVolume(bdys0, 0, lsc, IoR,
                                                       abs_matrix,
                                                       wave_len_min_matrix,
                                                       wave_len_max_matrix))
            # add bottom surface
            lsc[0].bdy_list.append(lsccls.TransparentBoundary(bdys0[0], lsc[0]))

            # add right interface with film
            lsc[0].bdy_list.append(lsccls.TransparentBoundary(bdys0[1],
                                                              lsc[0]))
            # add interface with rest of matrix
            lsc[0].bdy_list.append(lsccls.TransparentBoundary(bdys0[2],
                                                              lsc[0]))
            # add left solar cell
            lsc[0].bdy_list.append(lsccls.PVBoundary(bdys0[3], lsc[0],
                                                     'diffuse', EQE_pv,
                                                     0, wave_len_min_pv,
                                                     wave_len_max_pv))
            # add front mirror
            lsc[0].bdy_list.append(lsccls.OpaqueBoundary(bdys0[4], lsc[0],
                                                         'specular', .05))
            # add back mirror
            lsc[0].bdy_list.append(lsccls.OpaqueBoundary(bdys0[5], lsc[0],
                                                         'specular', .05))

            # define dimensions/characteristics of Volume 1 - phosphor film

            # input boundaries by setting boundary points of each
            bdy1a = [[0, L1, H1_0], [0, L2, H2_0], [W, L2, H2_0], [W, L1, H1_0]]
            bdy1b = [[0, L2, H2_0], [0, L2, H2_1], [W, L2, H2_1], [W, L2, H2_0]]
            bdy1c = [[0, L1, H1_1], [0, L2, H2_1], [W, L2, H2_1], [W, L1, H1_1]]
            bdy1d = [[0, L1, H1_0], [0, L1, H1_1], [W, L1, H1_1], [W, L1, H1_0]]
            bdy1e = [[0, L1, H1_0], [0, L2, H2_0], [0, L2, H2_1], [0, L1, H1_1]]
            bdy1f = [[W, L1, H1_0], [W, L2, H2_0], [W, L2, H2_1], [W, L1, H1_1]]
            bdys1 = [bdy1a, bdy1b, bdy1c, bdy1d, bdy1e, bdy1f]

            lsc.vol_list.append(lsccls.ParticleVolume(bdys1, 1, lsc, IoR,
                                                      abs_matrix, particle,
                                                      wave_len_min_matrix,
                                                      wave_len_max_matrix))
            # add bottom surface
            lsc[1].bdy_list.append(lsccls.TransparentBoundary(bdys1[0], lsc[1]))

            # add right mirror
            lsc[1].bdy_list.append(lsccls.OpaqueBoundary(bdys1[1], lsc[1],
                                                         'specular', .05))
            # add top surface
            lsc[1].bdy_list.append(lsccls.TransparentBoundary(bdys1[2], lsc[1]))
            # add left interface with mirror gap
            lsc[1].bdy_list.append(lsccls.TransparentBoundary(bdys1[3], lsc[1]))
            # add front mirror
            lsc[1].bdy_list.append(lsccls.OpaqueBoundary(bdys1[4], lsc[1],
                                                         'specular', .05))
            # add back mirror
            lsc[1].bdy_list.append(lsccls.OpaqueBoundary(bdys1[5], lsc[1],
                                                         'specular', .05))

            # define dimensions/characteristics of Volume 2 - rest of matrix

            # input boundaries by setting boundary points of each
            bdy2a = [[0, L0, H0_1], [0, L2, H2_1], [W, L2, H2_1], [W, L0, H0_1]]
            bdy2b = [[0, L2, H2_1], [0, L2, H2_2], [W, L2, H2_2], [W, L2, H2_1]]
            bdy2c = [[0, L0, H0_2], [0, L2, H2_2], [W, L2, H2_2], [W, L0, H0_2]]
            bdy2d = [[0, L0, H0_1], [0, L0, H0_2], [W, L0, H0_2], [W, L0, H0_1]]
            bdy2e = [[0, L0, H0_1], [0, L2, H2_1], [0, L2, H2_2], [0, L0, H0_2]]
            bdy2f = [[W, L0, H0_1], [W, L2, H2_1], [W, L2, H2_2], [W, L0, H0_2]]
            bdys2 = [bdy2a, bdy2b, bdy2c, bdy2d, bdy2e, bdy2f]

            # define volume
            lsc.vol_list.append(lsccls.AbsorbingVolume(bdys2, 2, lsc, IoR,
                                                       abs_matrix,
                                                       wave_len_min_matrix,
                                                       wave_len_max_matrix))
            # add interface with mirror gap and phosphor film
            lsc[2].bdy_list.append(lsccls.TransparentBoundary(bdys2[0],
                                                              lsc[2]))
            # add right mirror
            lsc[2].bdy_list.append(lsccls.OpaqueBoundary(bdys2[1], lsc[2],
                                                         'specular', .05))
            # add top surface
            lsc[2].bdy_list.append(lsccls.TransparentBoundary(bdys2[2],
                                                              lsc[2]))
            # add solar cell
            lsc[2].bdy_list.append(
                lsccls.PVBoundary(bdys2[3], lsc[2], 'diffuse', EQE_pv, 0,
                                  wave_len_min_pv, wave_len_max_pv))
            # add front mirror
            lsc[2].bdy_list.append(lsccls.OpaqueBoundary(bdys2[4], lsc[2],
                                                         'specular', .05))
            # add back mirror
            lsc[2].bdy_list.append(lsccls.OpaqueBoundary(bdys2[5], lsc[2],
                                                         'specular', .05))

            # Prepare data inputs for LSC simulation
            if light_form == 'direct' and results == 'matrix':
                theta_o = df.index[i]
                phi_o = df.columns[j]

            if light_form != 'direct' and results == 'matrix':
                tilt = df.index[i]

            lsc.matching_pairs()
            I = Einc * math.cos(math.radians(theta_o)) * (L * W)
            theta_o = math.radians(theta_o + 180)  # adjust theta to head down
            phi_o = math.radians(phi_o + 90)  # adjust phi
            tilt = math.radians(tilt)  # adjust tilt

            # Run LSC trials, determining fate of every bundle
            starting_vol = len(lsc) - 1
            starting_bdy = 2

            lsc.main(trials, L, W, H, light_form, theta_o,
                     phi_o, tilt, starting_vol, starting_bdy, I,
                     emi_source, emi_source_max, particle)

            # Process data outputs from all LSC trials

            # determine if all bundles in volume 0 are accounted for
            errorcount0 = (lsc[0].bundles_absorbed +
                           lsc[0][0].bundles_absorbed +
                           lsc[0][1].bundles_reflected +
                           lsc[0][1].bundles_refracted +
                           lsc[0][2].bundles_reflected +
                           lsc[0][2].bundles_refracted +
                           lsc[0][3].bundles_absorbed +
                           lsc[0][4].bundles_absorbed +
                           lsc[0][5].bundles_absorbed)

            # determine if all bundles in volume 1 are accounted for
            errorcount1 = (lsc[1].bundles_absorbed +
                           lsc[1][0].bundles_absorbed +
                           lsc[1][1].bundles_absorbed +
                           lsc[1][2].bundles_reflected +
                           lsc[1][2].bundles_refracted +
                           lsc[1][3].bundles_reflected +
                           lsc[1][3].bundles_refracted +
                           lsc[1][4].bundles_absorbed +
                           lsc[1][5].bundles_absorbed +
                           particle.bundles_absorbed)

            # determine if all bundles in volume 2 are accounted for
            errorcount2 = (lsc[2].bundles_absorbed +
                           lsc[2][0].bundles_reflected +
                           lsc[2][0].bundles_refracted +
                           lsc[2][1].bundles_absorbed +
                           lsc[2][2].bundles_reflected +
                           lsc[2][2].bundles_refracted +
                           lsc[2][3].bundles_absorbed +
                           lsc[2][4].bundles_absorbed +
                           lsc[2][5].bundles_absorbed)

            error = (errorcount0 + errorcount1 + errorcount2) / trials

            if error != 1:
                print("\nENERGY IS NOT CONSERVED!!!!!")

            if light_form == 'direct' and results == 'matrix':
                df.iloc[i, j] = lsc
                df_oe.iloc[i, j] = lsc.optical_efficiency
                df_m.iloc[i, j] = lsc.m

            if light_form != 'direct' and results == 'matrix':
                df.iloc[i, 0] = lsc
                df_oe.iloc[i, 0] = lsc.optical_efficiency
                df_m.iloc[i, 0] = lsc.m

    if results == 'matrix':
        df_params.to_csv('mc_params.csv')
        lsc_outputs = {'parameters': mc_params, 'opt_eff': df_oe,
                       'spect_mismatch': df_m, 'object': df, }
        if light_form == 'direct':
            df_oe.to_csv('opt_eff_beam_circ.csv')
            df_m.to_csv('spect_mismatch_beam_circ.csv')
        if light_form == 'diffuse':
            df_oe.to_csv('opt_eff_iso.csv')
            df_m.to_csv('spect_mismatch_iso.csv')
        if light_form == 'ground':
            df_oe.to_csv('opt_eff_grnd.csv')
            df_m.to_csv('spect_mismatch_grnd.csv')

    else:
        print(time.time() - start_time)
        print(light_form)
        print(lsc.optical_efficiency)
        print(lsc.m)
        print(tilt)
        lsc_outputs = {'parameters': mc_params, 'opt_eff': lsc.optical_efficiency,
                       'spect_mismatch': lsc.m, 'object': lsc}

    return lsc_outputs

#if __name__ == "__main__":
    #lsc_outputs=wedge(100,theta_o=60,phi_o=0.001)
    #lsc_outputs = wedge(100, theta_o=0.001, phi_o=0.001)
    #lsc_outputs = wedge(10000, phi_o=0.001, theta_o=180)
    #print(lsc_outputs)