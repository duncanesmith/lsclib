import os
import math
import random
import time
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import PolynomialFeatures 
from sklearn.linear_model import LinearRegression
from scipy.interpolate import InterpolatedUnivariateSpline



def regression(directory):
    f = []
    for filename in os.listdir(directory):
        if filename.endswith(".xlsx"): 
            f.append(filename)
            continue
        else:
            continue
    theta_equations = [] 
    phi_equations = []
    for j in range(0,len(f)):
        df = pd.read_excel(f[j])
        vectx = df["VectorX"]
        vecty = df['VectorY']
        vectz = df['VectorZ']
       
        vectx = list(vectx)
        vecty = list(vecty)
        vectz = list(vectz)
        
        theta = []
        phi = []
        for i in range(0,len(vectx)):
            p = math.atan(vecty[i]/vectx[i]) 
            phi.append(p)
            t = math.atan(math.sqrt(vectx[i]**2 + vecty[i]**2)/vectz[i])+math.pi
            theta.append(t)
        index = np.linspace(1,100,100)
        spl2 = InterpolatedUnivariateSpline(index,phi)
        spl = InterpolatedUnivariateSpline(index,theta) #equation
        theta_equations.append(spl)
        phi_equations.append(spl2)
        #multivariate linear regression
    
    d = pd.DataFrame(columns = np.linspace(0,85,18), index = np.linspace(400,1200,17))
    p = pd.DataFrame(columns = np.linspace(0,85,18), index = np.linspace(400,1200,17))
    p2 = pd.DataFrame(columns = np.linspace(0,85,18), index = np.linspace(400,1200,17))
    d2 = pd.DataFrame(columns = np.linspace(0,85,18), index = np.linspace(400,1200,17))
    k = 0
    for i in range(0,len(d)):
        for j in range(0,18):
            d.iloc[i,j] = theta_equations[k]
            p.iloc[i,j] = phi_equations[k]
           #d2.iloc[i,j] = theta_equations[k].__call__([np.linspace(1,100,10)])
            #p2.iloc[i,j] = phi_equations[k].__call__([np.linspace(1,100,10)])
            k+=1
    return [d,p]
files = regression(r"C:\Users\vinny\Downloads\Regression Excel Sheets (3)\Regression Excel Sheets")


def equations(wavelength,angle,theta_df, phi_df): 
    p = phi_df
    d = theta_df
    if wavelength%100 != 0 or wavelength%100 != 50: # round to the nearest 50
        wavelength = int(50 * round(float(wavelength)/50))
    
    if angle%10 !=0 or angle%5 != 5:
        angle = int(5*round(float(angle)/5))
    
    rand = random.uniform(1,100)
    #phi = p.loc[wavelength,angle].__call__(rand)
    #theta = d.loc[wavelength,angle].__call__(rand)
    phi = p.loc[wavelength,angle]
    theta = d.loc[wavelength,angle]
    
    #theta_linreg = model.predict([rand,wavelength,angle]) #Add model to code
    #xs = np.linspace(1,100,100)
    return [theta,phi] #one function

#start = time.time()
wavelength = 900
theta_input = 70
trial = equations(500, 50, files[0],files[1])
#end = time.time()
#print(end-start)

theta_eq = (trial[0])
phi_eq = (trial[1])
a = []
theta_linreg = []
phi_linreg = []
for i in range(0,1000):
    b = (random.uniform(1,100))
    a.append(b)



fig,ax= plt.subplots(2,1, sharex=False, sharey=False)
ax[0].hist(theta_eq(a), bins=20)
#ax[0].set(xlabel="index")
ax[0].set(ylabel="Theta")
ax[1].hist(phi_eq(a), bins=20)
#ax[1].set(xlabel="index")
ax[1].set(ylabel="phi")

#directory = r"C:\Users\patelb6\Downloads\Senior Year\Research\Regression Excel Sheets"
#f = []
#for filename in os.listdir(directory):
    #if filename.endswith(".xlsx"): 
        #f.append(filename)
        #continue
    #else:
        #continue
#theta_equations = [] 
#phi_equations = []
#for j in range(0,len(f)):
    #df = pd.read_excel(f[j])
    #vectx = df["VectorX"]
    #vecty = df['VectorY']
    #vectz = df['VectorZ']
   
    #vectx = list(vectx)
    #vecty = list(vecty)
    #vectz = list(vectz)
    
    #theta = []
    #phi = []
    #for i in range(0,len(vectx)):
        #p = math.atan(vecty[i]/vectx[i]) 
        #phi.append(p)
        #t = math.atan(math.sqrt(vectx[i]**2 + vecty[i]**2)/vectz[i])
        #theta.append(t)
    #index = np.linspace(1,100,100)
    #spl2 = InterpolatedUnivariateSpline(index,phi)
    #spl = InterpolatedUnivariateSpline(index,theta) #equation
    #theta_equations.append(spl)
    #phi_equations.append(spl2)
    ##multivariate linear regression

#d = pd.DataFrame(columns = np.linspace(0,90,19), index = np.linspace(400,1200,17))
#p = pd.DataFrame(columns = np.linspace(0,90,19), index = np.linspace(400,1200,17))
#p2 = pd.DataFrame(columns = np.linspace(0,90,19), index = np.linspace(400,1200,17))
#d2 = pd.DataFrame(columns = np.linspace(0,90,19), index = np.linspace(400,1200,17))
#k = 0
#for i in range(0,len(d)):
    #for j in range(0,19):
        #d.iloc[i,j] = theta_equations[k]
        #p.iloc[i,j] = phi_equations[k]
        #d2.iloc[i,j] = theta_equations[k].__call__([np.linspace(1,100,10)])
        #p2.iloc[i,j] = phi_equations[k].__call__([np.linspace(1,100,10)])
        #k+=1



#wavelength = 600
#angle = 45
#if wavelength%100 != 0 or wavelength%100 != 50: # round to the nearest 50
    #wavelength = int(50 * round(float(wavelength)/50))

#if angle%10 !=0 or angle%5 != 5: #rounds to the nearest 10
    #angle = int(5*round(float(angle)/5))



#phi = p.loc[wavelength,angle]
#theta = d.loc[wavelength,angle]



#a = []
#theta_linreg = []
#phi_linreg = []
#for i in range(0,1000):
    #b = (random.uniform(1,100))
    #a.append(b)


##fig,ax= plt.subplots(2,2, sharex=False, sharey=False)
##ax[0,0].hist(theta_eq(a), bins=20)
###ax[0].set(xlabel="index")
##ax[0,0].set(ylabel="Theta")
##ax[0,1].hist(phi_eq(a), bins=20)
###ax[1].set(xlabel="index")
##ax[0,1].set(ylabel="phi")