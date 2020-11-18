# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 11:15:45 2019

@author: smithd24
"""
import pandas as pd
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline


def excel_read():

    # read silicone matrix absorption spectrum data from excel
    abs_matrix = pd.read_excel('absorption_matrix.xlsx')
    # read pv external quantum efficiency spectrum data from excel
    EQE_pv = pd.read_excel('eqe_pv.xlsx')
    # read pv internal quantum efficiency spectrum data from excel
    IQE_pv = pd.read_excel('iqe_pv.xlsx')
    # read xenon emission spectrum data from excel
    emi_source = pd.read_excel('emission_source.xlsx')
    # read phosphor absorption spectrum data from excel
    abs_particle = pd.read_excel('absorption_particle.xlsx')    
    # read phosphor emission spectrum data from excel
    emi_particle = pd.read_excel('emission_particle.xlsx')
    
    return abs_matrix, EQE_pv, IQE_pv, emi_source, abs_particle, emi_particle

def spline(dataset):
  "Creates a polynomial regression model for the given degree"
  
  num_cols = np.size(dataset,1)
  X = dataset.iloc[:, 0:1].values
  Y = dataset.iloc[:, num_cols - 1].values
  max_value_x = max(X)[0]
  spl = InterpolatedUnivariateSpline(X, Y)
  
  return spl, max_value_x