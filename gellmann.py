#!/usr/bin/env python
# coding: utf-8
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 12:02:08 2019

@author: Rohan
"""
import numpy as np
from sympy import *
from IPython.display import display
init_printing(use_latex=True,forecolor="White",fontsize="10pt")

lam=np.zeros((9,3,3)).astype(complex)

lam[1,0,1],lam[1,1,0]=1,1
lam[2,0,1],lam[2,1,0]=-1j,1j
lam[3,0,0],lam[3,1,1]=1,-1
lam[4,0,2],lam[4,2,0]=1,1
lam[5,0,2],lam[5,2,0]=-1j,1j
lam[6,1,2],lam[6,2,1]=1,1
lam[7,1,2],lam[7,2,1]=-1j,1j
lam[8,0,0],lam[8,1,1],lam[8,2,2]=1/sqrt(3),1/sqrt(3),-2/sqrt(3)
            
R=np.zeros((9,9,3,3)).astype(complex)
for i in range(8):
    for j in range(8):
        R[i+1,j+1,:,:]=Matrix(np.matmul(lam[i+1,:,:],lam[j+1,:,:]))
    
def return_lamda():
    global lam
    return lam

def struc_cons():
    global lam
    f=np.zeros((9,9,9)).astype(complex)
    d=np.zeros((9,9,9)).astype(complex)
    for i in range(8):
        for k in range(8):
            for l in range(8):
                f[i+1,k+1,l+1]=(-1j/4)*np.trace(np.matmul(lam[i+1,:,:],(np.matmul(lam[k+1,:,:],lam[l+1,:,:])-np.matmul(lam[l+1,:,:],lam[k+1,:,:]))))
                d[i+1,k+1,l+1]=(1/4)*np.trace(np.matmul(lam[i+1,:,:],(np.matmul(lam[k+1,:,:],lam[l+1,:,:])+np.matmul(lam[l+1,:,:],lam[k+1,:,:]))))
    return np.real(f)

if __name__ != "__main__":
    struc_cons()

