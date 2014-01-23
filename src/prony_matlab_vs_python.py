#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Jan 22 16:56 2014

@author: Sammy Pfeiffer
@email: sammypfeiffer@gmail.com

Test our python prony vs the MatlabTM prony
"""
from matlabpipe import MatlabPipe
import numpy as np
from prony import prony

OKGREEN = '\033[92m'
FAIL = '\033[91m'
ENDC = '\033[0m'

if __name__ == '__main__':
    matlab = MatlabPipe(matlab_version='2013a')
    matlab.open()
    
    cmd = """[b,a] = butter(4,0.2);
    impulseResp = impz(b,a);
    denOrder=4;
    numOrder=4;
    [Num,Den]=prony(impulseResp,numOrder,denOrder);"""
    out = matlab.eval(cmd)
    # Get matlab inputs
    matlab_impulseResp = matlab.get('impulseResp')
    matlab_denOrder = matlab.get('denOrder')
    matlab_numOrder = matlab.get('numOrder')
    # Get matlab output of prony
    matlab_Num = matlab.get('Num')
    matlab_Den = matlab.get('Den')
     
    # Use same inputs on python prony
    [py_Num, py_Den] = prony(matlab_impulseResp,matlab_numOrder,matlab_denOrder)
     
    # Check if we get the same result
    print "\nMatlab Num vs Python Num"
    if np.allclose(matlab_Num, py_Num): # see if they are equal with relative tolerance=1e-05 and abs tol=1e-08
        print OKGREEN + "Coincident results!" + ENDC
    else:
        print FAIL +"Different results..." + ENDC
    print matlab_Num
    print py_Num
     
    print "\nMatlab Den vs Python Den"
    if np.allclose(matlab_Den, py_Den):
        print OKGREEN + "Coincident results!" + ENDC
    else:
        print FAIL +"Different results..." + ENDC
    print matlab_Den
    print py_Den
    
    
    
    