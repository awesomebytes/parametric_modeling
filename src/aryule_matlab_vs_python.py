#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Jan 22 22:26 2014

@author: Sammy Pfeiffer
@email: sammypfeiffer@gmail.com

Test our python aryule vs the MatlabTM aryule

MATLAB testing code from:
http://www.music.mcgill.ca/~gary/307/matlab/predict.m
"""
from matlabpipe import MatlabPipe
import numpy as np
from aryule import aryule

OKGREEN = '\033[92m'
FAIL = '\033[91m'
ENDC = '\033[0m'

if __name__ == '__main__':
    matlab = MatlabPipe(matlab_version='2013a')
    matlab.open()
    
    cmd = """A=[1 -2.7607 3.8106 -2.6535 0.9238];
        y=filter(1,A,0.2*randn(1024,1));
        [ar_coeffs, e, k]=aryule(y,4);"""
    out = matlab.eval(cmd)
    # Get matlab inputs
    matlab_A = matlab.get('A')
    matlab_y = matlab.get('y')
    # Get matlab output of aryule
    matlab_ar_coeffs = matlab.get('ar_coeffs')
    matlab_e = matlab.get('e')
    matlab_k = matlab.get('k')

     
    # Use same inputs on python aryule
    [py_ar_coeffs, py_e, py_k] = aryule(matlab_y, 4)
     
    # Check if we get the same result
    print "\nMatlab ar_coeffs vs Python ar_coeffs"
    if np.allclose(matlab_ar_coeffs, py_ar_coeffs): # see if they are equal with relative tolerance=1e-05 and abs tol=1e-08
        print OKGREEN + "Coincident results!" + ENDC
    else:
        print FAIL +"Different results..." + ENDC
    print matlab_ar_coeffs
    print py_ar_coeffs
     
    
    print "\nMatlab e vs Python e"
    if np.allclose(matlab_e, py_e): # see if they are equal with relative tolerance=1e-05 and abs tol=1e-08
        print OKGREEN + "Coincident results!" + ENDC
    else:
        print FAIL +"Different results..." + ENDC
    print matlab_e
    print py_e
    
    print "\nMatlab k vs Python k"
    if np.allclose(matlab_k, py_k): # see if they are equal with relative tolerance=1e-05 and abs tol=1e-08
        print OKGREEN + "Coincident results!" + ENDC
    else:
        print FAIL +"Different results..." + ENDC
    print matlab_k
    print py_k
    