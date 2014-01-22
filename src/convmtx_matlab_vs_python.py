#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Jan 22 17:18 2014

@author: Sammy Pfeiffer

Test our python cnvmtx vs the MatlabTM cnvmtx
"""
from matlabpipe import MatlabPipe
import numpy as np
from convmtx import convmtx

OKGREEN = '\033[92m'
FAIL = '\033[91m'
ENDC = '\033[0m'

if __name__ == '__main__':
    matlab = MatlabPipe(matlab_version='2013a')
    matlab.open()
    
    cmd = """h = [1 2 3 2 1];
    c_matrix = convmtx(h,7);"""
    out = matlab.eval(cmd)
    # Get matlab inputs
    matlab_h = matlab.get('h')
    # Get matlab output of convmtx
    matlab_c_matrix = matlab.get('c_matrix')

     
    # Use same inputs on python convmtx
    vertical_h = np.vstack((matlab_h))
    py_c_matrix = convmtx(vertical_h,7)
     
    # Check if we get the same result
    print "\nMatlab c_matrix vs Python c_matrix"
    if np.allclose(matlab_c_matrix, py_c_matrix): # see if they are equal with relative tolerance=1e-05 and abs tol=1e-08
        print OKGREEN + "Coincident results!" + ENDC
    else:
        print FAIL +"Different results..." + ENDC
    print matlab_c_matrix
    print py_c_matrix
     
    
    