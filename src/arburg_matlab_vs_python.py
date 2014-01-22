#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Jan 22 17:18 2014

@author: Sammy Pfeiffer

Test our python arburg vs the MatlabTM arburg
"""
from matlabpipe import MatlabPipe
import numpy as np
from arburg import arburg

OKGREEN = '\033[92m'
FAIL = '\033[91m'
ENDC = '\033[0m'

if __name__ == '__main__':
    matlab = MatlabPipe(matlab_version='2013a')
    matlab.open()
    
    cmd = """input = [0 1 2 3 4 5]
    [ar_coeffs, NoiseVariance] = arburg(input, 4)"""
    out = matlab.eval(cmd)
    # Get matlab inputs
    matlab_input = matlab.get('input')
    # Get matlab output of stmcb
    matlab_ar_coeffs = matlab.get('ar_coeffs')
    matlab_NoiseVariance = matlab.get('NoiseVariance')
    #matlab_sys_model = matlab.get('sys_model')
     
    # Use same inputs on python stmcb
    [py_ar_coeffs, py_NoiseVariance] = arburg(matlab_input, 4)
     
    # Check if we get the same result
    print "\nMatlab ar_coeffs vs Python ar_coeffs"
    if np.allclose(matlab_ar_coeffs, py_ar_coeffs): # see if they are equal with relative tolerance=1e-05 and abs tol=1e-08
        print OKGREEN + "Coincident results!" + ENDC
    else:
        print FAIL +"Different results..." + ENDC
    print matlab_ar_coeffs
    print py_ar_coeffs
     
    print "\nMatlab NoiseVariance vs Python NoiseVariance"
    if np.allclose(matlab_NoiseVariance, py_NoiseVariance):
        print OKGREEN + "Coincident results!" + ENDC
    else:
        print FAIL +"Different results..." + ENDC
    print matlab_NoiseVariance
    print py_NoiseVariance
    
    
    
    