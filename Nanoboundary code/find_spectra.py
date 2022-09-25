# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 09:40:28 2020

@author: smithd24
"""

import numpy as np
import lsc_main as lm
import time
import pandas as pd

theta = [0.000001]
flux_gain = np.zeros(len(theta))
simulations = np.ones(len(theta))*2500000
total_time = 0

for i in range(0,len(flux_gain)):
    
    start_time = time.time()
    simulations_int = int(simulations[i])
    print('\n')
    print(simulations_int)
    print(theta[i])
    LSC = lm.LSC_wedge_main(simulations_int, theta_o = theta[i])
    elapsed_time = time.time() - start_time
    
    wave_len_log = pd.DataFrame(LSC.wave_len_log)
    norm_inc_spectrum = LSC.norm_inc_spectrum
    wave_len_log_pv = pd.DataFrame(LSC.wave_len_log_pv)
    norm_inc_spectrum_pv = LSC.norm_inc_spectrum_pv
    
    flux_gain[i] = (LSC.Isc_cell/(0.022*0.0045))/(LSC.Isc_i/(0.05*0.022)) 
    total_time += elapsed_time
    print(flux_gain[i])
    print(i)
    print(elapsed_time)
    print(total_time)

spectra_binned = pd.concat([norm_inc_spectrum, norm_inc_spectrum_pv], axis=1)
spectra_binned.plot()
writer = pd.ExcelWriter('spectra_5bundles4.5height.xlsx', engine = 'xlsxwriter')
spectra_binned.to_excel(writer,'Sheet1')
wave_len_log_pv.to_excel(writer,'Sheet2')
writer.save()