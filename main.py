# -*- coding: utf-8 -*-
"""
Created on Mon May 16 11:03:37 2022

@author: D krishnathejus
"""

import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np

def enclosure(Ra,Pr,Length,Grid):
    theta = 1-1/Grid*np.meshgrid(np.arange(0,Grid),np.arange(0,Grid))[0]
    psi=np.zeros([Grid,Grid])
    omega=np.zeros([Grid,Grid])
    dx=dy=Length/Grid
    for row in range (Grid):
        theta[row,0]=1
        theta[row,Grid-1]=0
    Nu_ca = 0
    Nu_ha = 1
    iter = 0
    
    while(abs(Nu_ca - Nu_ha)>0.000001 and iter < 1000):
           
    
        for i in range(1,Grid-1):
            for j in range(1,Grid-1):
                psi[i,j] += 0.5*((-omega[i,j]-(psi[i+1,j]+psi[i-1,j])/(dy*dy)-(psi[i,j+1]+psi[i,j-1])/(dx*dx))/((-2/(dx*dx))+(-2/(dy*dy))) -  psi[i,j])
                omega[i,j] += 0.5*((((psi[i+1,j]-psi[i-1,j])/(2*dy))*(omega[i,j+1]-omega[i,j-1])/(2*dx) - ((psi[i,j+1]-psi[i,j-1])/(2*dx))*(omega[i+1,j]-omega[i-1,j])/(2*dy) - Ra*Pr*(theta[i,j+1]-theta[i,j-1])/(2*dx) - Pr*((omega[i,j+1]+omega[i,j-1])/(dx*dx)+(omega[i+1,j]+omega[i-1,j])/(dy*dy)))/(Pr*((-2/(dx*dx))+(-2/(dy*dy)))) -  omega[i,j])
                theta[i,j] += 0.5*((((psi[i+1,j]-psi[i-1,j])/(2*dy))*(theta[i,j+1]-theta[i,j-1])/(2*dx) - ((psi[i,j+1]-psi[i,j-1])/(2*dx))*(theta[i+1,j]-theta[i-1,j])/(2*dy) - (theta[i,j+1]+theta[i,j-1])/(dx*dx) - (theta[i+1,j]+theta[i-1,j])/(dy*dy))/((-2/(dx*dx))+(-2/(dy*dy))) -  theta[i,j])
      
        for i in range (Grid):
            omega[i,0] = -2*(psi[i,1])/(dx*dx)
            omega[i,Grid-1]= -2*(psi[i,Grid-2])/(dx*dx)
    
        for i in range(1,Grid-1):
            omega[0,i]= -2*(psi[1,i])/(dy*dy)
            omega[Grid-1,i]= -2*(psi[Grid-2,i])/(dy*dy)
            
            theta[0,i]=(4/3)*(theta[1,i] - theta[2,i]/4)
            theta[Grid-1,i]=(4/3)*(theta[Grid-2,i] - theta[Grid-3,i]/4)
    
        
        #Calculating Nusselts Number
        Nu_h = [(1- theta[row,1])/dx for row in range(Grid)]
        Nu_c = [(theta[row,Grid-2])/dx for row in range(Grid)]
        Nu_ha = integrate.simps(Nu_h,dx = dx)
        Nu_ca = integrate.simps(Nu_c,dx = dx)
        


        
        iter += 1
    
    
    #streamline plot
    fig = plt.figure(figsize=(15,12))
    fig, ax = plt.subplots(figsize=(15,12))
    ax.tick_params(axis='both', which='major', labelsize=20)
    plt.title("Stream Function",fontdict={'fontsize':25})
    fig.suptitle("Ra="+str(Ra))
    Streamfunction = plt.contourf(psi)
    plt.clabel(Streamfunction, inline=0, fontsize=0)
    plt.set_cmap('plasma')
    size_label= fig.colorbar(Streamfunction, ax=ax)
    size_label.ax.tick_params(labelsize=20)
    plt.axis('off')
    plt.show()
    
    #Vorticity Field
    fig = plt.figure(figsize=(15,12))
    fig, ax = plt.subplots(figsize=(15,12))
    ax.tick_params(axis='both', which='major', labelsize=20)
    fig.suptitle("Ra="+str(Ra))
    Vorticity = plt.contourf(omega)
    plt.clabel(Vorticity, inline=0, fontsize=0)
    plt.set_cmap('ocean')
    size_label= fig.colorbar(Vorticity, ax=ax)
    size_label.ax.tick_params(labelsize=20)
    plt.axis('off')
    plt.show()

    #Isotherm
    fig = plt.figure(figsize=(15,12))
    fig, ax = plt.subplots(figsize=(15,12))
    ax.tick_params(axis='both', which='major', labelsize=20)
    fig.suptitle("Ra="+str(Ra))
    plt.xlabel('Isotherm')
    plt.ylabel('Y-Axis')
    Isotherm = plt.contourf(theta)
    plt.clabel(Isotherm, inline=0, fontsize=0)
    plt.set_cmap('autumn')
    size_label= fig.colorbar(Isotherm, ax=ax)
    size_label.ax.tick_params(labelsize=20)
    plt.axis('off')
    plt.show()
    
    
enclosure(1000,0.71,1,110) #function call for Rayleigh number 1000
enclosure(10000,0.71,1,110) #function call for Rayleigh number 10000
enclosure(100000,0.71,1,110) #function call for Rayleigh number 100000
