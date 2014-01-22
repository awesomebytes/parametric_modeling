#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Jan 22 17:18 2014

@author: Sammy Pfeiffer

Test our python stmcb vs the MatlabTM stmcb
"""
from matlabpipe import MatlabPipe
import numpy as np
from stmcb import stmcb

OKGREEN = '\033[92m'
FAIL = '\033[91m'
ENDC = '\033[0m'

if __name__ == '__main__':
    matlab = MatlabPipe(matlab_version='2013a')
    matlab.open()
    
    cmd = """[b,a] = butter(6,0.2);
    h = filter(b,a,[1 zeros(1,100)]);
    [bb,aa] = stmcb(h,4,4);"""
    out = matlab.eval(cmd)
    # Get matlab inputs
    matlab_b = matlab.get('b')
    matlab_a = matlab.get('a')
    matlab_h = matlab.get('h')
    # Get matlab output of stmcb
    matlab_bb = matlab.get('bb')
    matlab_aa = matlab.get('aa')
     
    # Use same inputs on python stmcb
    [py_bb, py_aa] = stmcb(matlab_h,4,4)
     
    # Check if we get the same result
    print "\nMatlab bb vs Python bb"
    if np.allclose(matlab_bb, py_bb): # see if they are equal with relative tolerance=1e-05 and abs tol=1e-08
        print OKGREEN + "Coincident results!" + ENDC
    else:
        print FAIL +"Different results..." + ENDC
    print matlab_bb
    print py_bb
     
    print "\nMatlab matlab_aa vs Python py_aa"
    if np.allclose(matlab_aa, py_aa):
        print OKGREEN + "Coincident results!" + ENDC
    else:
        print FAIL +"Different results..." + ENDC
    print matlab_aa
    print py_aa
    
    
    
    