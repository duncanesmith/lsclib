# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 11:15:45 2019

@author: smithd24
"""
import pandas as pd
import numpy as np
import operator
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from sklearn.metrics import mean_squared_error, r2_score


def LSC_excel_read():

    # read silicone matrix absorption spectrum data from excel
    abs_silicone = pd.read_excel('absorption_silicone.xlsx')
    abs_silicone = abs_silicone.values  # convert dataframe to array using numpy
    
    # read pv absorption spectrum data from excel
    abs_pv = pd.read_excel('absorption_pv.xlsx')
    abs_pv = abs_pv.values  # convert dataframe to array using numpy

    # read xenon emission spectrum data from excel
    emi_xenon = pd.read_excel('emission_xenon.xlsx')
    emi_xenon = emi_xenon.values  # convert dataframe to array using numpy
    
    # read phosphor absorption spectrum data from excel
    abs_phosphor = pd.read_excel('absorption_phosphor.xlsx')
    abs_phosphor = abs_phosphor.values
    
    # read phosphor emission spectrum data from excel
    emi_phosphor = pd.read_excel('emission_phosphor.xlsx')
    emi_phosphor = emi_phosphor.values # convert dataframe to array using numpy
    
    return abs_silicone, abs_pv, emi_xenon, abs_phosphor, emi_phosphor

def LSC_excel_read_spl():

    # read silicone matrix absorption spectrum data from excel
    abs_silicone = pd.read_excel('absorption_silicone.xlsx')
    # read pv external quantum efficiency spectrum data from excel
    EQE_pv = pd.read_excel('eqe_pv.xlsx')
    # read pv internal quantum efficiency spectrum data from excel
    IQE_pv = pd.read_excel('iqe_pv.xlsx')
    # read xenon emission spectrum data from excel
    emi_xenon = pd.read_excel('emission_xenon.xlsx')
    # read phosphor absorption spectrum data from excel
    abs_phosphor = pd.read_excel('absorption_phosphor.xlsx')    
    # read phosphor emission spectrum data from excel
    emi_phosphor = pd.read_excel('emission_phosphor.xlsx')
    
    return abs_silicone, EQE_pv, IQE_pv, emi_xenon, abs_phosphor, emi_phosphor

def spline(dataset, smoothing):
  "Creates a polynomial regression model for the given degree"
  
  num_cols = np.size(dataset,1)
  X = dataset.iloc[:, 0:1].values
  Y = dataset.iloc[:, num_cols - 1].values
  max_value_x = max(X)[0]
  
  spl = UnivariateSpline(X, Y)
  spl.set_smoothing_factor(smoothing)
  
  # Y_predict = spl.__call__(X)
  
  # # evaluating the model on training dataset
  # rmse = np.sqrt(mean_squared_error(Y, Y_predict))
  # r2 = r2_score(Y, Y_predict)
  
  
  # print("The model performance for the training set")
  # print("-------------------------------------------")
  # print("RMSE of training set is {}".format(rmse))
  # print("R2 score of training set is {}".format(r2))
      
      
  # plt.scatter(X, Y, s=10)
  # sort_axis = operator.itemgetter(0)
  # sorted_zip = sorted(zip(X,Y_predict), key=sort_axis)
  # X, y_predict = zip(*sorted_zip)
  # plt.plot(X, Y_predict, color='m')
  # plt.show()
  # print(spl.get_knots())
  
  return spl, max_value_x