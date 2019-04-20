#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 10:20:28 2019

@author: dtubin
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
kB=1.38064852*1e-23 #m^2 kg s^-2 K^-1
planck=6.62607004*1e-34 #m^2 kg  s^-1
c=299792458. #m  s^-1


def f(x, A, B): 
    return A*x + B


def fit(x,y):
    A,B = curve_fit(f, x, y)[0]
    diag=curve_fit(f, x, y)[1]
    A_err,B_err = np.sqrt(np.diag(diag))
    return A,A_err,B,B_err



def T_rot(file,emision,degeneracy,energy_u):
    deg=np.array(degeneracy) 
    freq=file[:,0]*1e9
    inte=file[:,1]*1000
    Nu=[]
    for i in range(0,len(freq)):
        Nu.append((8*np.pi*inte[i]*kB*(freq[i])**2)/(emision[i]*planck*(c*100)**3))
        #print (freq[i]/1e9,'&','{:.4e}'.format(Nu[i]))
    Nu=np.array(Nu)
    y_axis=(np.log(Nu/deg))
    x_axis=np.array(energy_u)
    slope,slope_err,interc,interc_err=fit(x_axis,y_axis)
    print(slope,slope_err,interc,interc_err)
    x=np.linspace(x_axis[0],x_axis[4],100)
    t_rot=np.abs(1/slope)
    print('{:.4e}'.format(t_rot), '{:.4e}'.format(inter))
    t_rot_err=np.abs(slope_err*1/(slope*slope))
    print('err',t_rot_err)
    f1=plt.figure()
    plt.plot(x_axis,y_axis, 'bo', label='Data')
    plt.plot(x,slope*x+interc, label=r'Fit: $log(\frac{N_{u}}{g_{u}})=A \cdot E_{u}+B$')
    plt.ylabel(r'$ln(\frac{N_{u}}{g_{u}})$',fontsize=20)
    plt.xlabel(r'$E_{u}[K]$', color='#1C2833',fontsize=20)
    plt.legend(loc='best')
    plt.title(r'Rotational Diagram')
    f1.savefig('hola.pdf')
    return t_rot,t_rot_err,interc,interc_err
    
def interpolate(temp,log_Z,t_rot):
    f=interp1d(temp,log_Z,kind='cubic',bounds_error=False,fill_value='extrapolate')

    x=np.linspace(0,300,1000)
    f2=plt.figure()
    ax2=f2.add_subplot(111)
    plt.plot(temp,log_Z,'ro',label='Table data')
    plt.plot(x,f(x),'g-',label=r'f(x): Interpolation')
    plt.ylabel(r'$log(Z_{H_{2}CO})$')
    plt.xlabel(r'Rotational Temperature[K]')
    plt.legend(loc='best')
    return f(t_rot)


file=np.loadtxt("/home/duuuuuuusan/Descargas/h2co-obs.dat")
deg=[15.,15.,21.,21.,27.]
Aul=10**np.array([-4.27564,-4.18926,-3.64396,-3.55754,-3.23040])
Eu=[21.92264, 22.61799, 32.05926, 33.44986,45.57057]

temp_rot,temp_rot_err,inter,inter_err=T_rot(file,Aul,deg,Eu)

temp_int=[300,225,150,75,37.5,18.75,9.375]
logZ=[3.4598,3.2724,3.0086,2.5584,2.1094,1.6501,1.1399]
partition_function=interpolate(temp_int,logZ,temp_rot)
N_h2c0=np.exp(inter+partition_function)
#N_h2co_err=np.sqrt(N_h2c0**2*inter_err**2)
print('{:.4e}'.format(N_h2c0))
#print(partition_function)
plt.show()