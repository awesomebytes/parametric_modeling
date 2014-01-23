#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Jan 22 16:56 2014

@author: Sammy Pfeiffer
@email: sammypfeiffer@gmail.com

Test our python invfreqs vs the MatlabTM invfreqs
"""
from matlabpipe import MatlabPipe
import numpy as np
from invfreqs import invfreqs

OKGREEN = '\033[92m'
FAIL = '\033[91m'
ENDC = '\033[0m'

if __name__ == '__main__':
    matlab = MatlabPipe(matlab_version='2013a')
    matlab.open()
    
    cmd = """b = [1 2 3 2 3];
    a = [1 2 3 2 1 4];
    [h,w] = freqs(b,a,64);
    [bb,aa] = invfreqs(h,w,4,5)"""
    out = matlab.eval(cmd)
    # Get matlab inputs
    matlab_b = matlab.get('b')
    matlab_a = matlab.get('a')
    matlab_h = matlab.get('h')
    matlab_w = matlab.get('w')
    # Get matlab output of invfreqs
    matlab_bb = matlab.get('bb')
    matlab_aa = matlab.get('aa')
     
    # Use same inputs on python invfreqs
    [py_bb, py_aa] = invfreqs(matlab_h, matlab_w, 4, 5)
     
    # Check if we get the same result
    print "\nMatlab bb vs Python bb"
    if np.allclose(matlab_bb, py_bb): # see if they are equal with relative tolerance=1e-05 and abs tol=1e-08
        print OKGREEN + "Coincident results!" + ENDC
    else:
        print FAIL +"Different results..." + ENDC
    print matlab_bb
    print py_bb
     
    print "\nMatlab Den vs Python Den"
    if np.allclose(matlab_aa, py_aa):
        print OKGREEN + "Coincident results!" + ENDC
    else:
        print FAIL +"Different results..." + ENDC
    print matlab_aa
    print py_aa
    
    
    
    