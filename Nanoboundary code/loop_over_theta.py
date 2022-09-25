# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 09:40:28 2020

@author: smithd24
"""

import numpy as np
import run
import time
import pandas as pd

theta = [.000001, .000001, .000001]
simulations = np.ones(len(theta))*800000
total_time = 0

LSC_stats = pd.DataFrame(columns = ['theta','optical_efficiency',
                                    'spectral_mismatch','Isc_i','Isc_cell'])

for i in range(0,len(theta)):
    
    start_time = time.time()
    simulations_int = int(simulations[i])
    print('\n')
    print(simulations_int)
    print(theta[i])
    LSC = run.wedge(simulations_int, Einc = 1000, theta_o = theta[i])
    elapsed_time = time.time() - start_time 
    new_row = {'theta': theta[i],
               'optical_efficiency': LSC.optical_efficiency,
               'spectral_mismatch': LSC.m,
               'Isc_i': LSC.Isc_i,
               'Isc_cell': LSC.Isc_cell}
    LSC_stats = LSC_stats.append(new_row, ignore_index = True)
    total_time += elapsed_time
    print(LSC.m)
    print(i)
    print(total_time)

stats = pd.DataFrame(LSC_stats)
writer = pd.ExcelWriter('stats800000planar3cm.xlsx')
stats.to_excel(writer,'Sheet1')
writer.save()