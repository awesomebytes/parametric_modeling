#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Jan 22 17:18 2014

@author: Sammy Pfeiffer
@email: sammypfeiffer@gmail.com

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
    
    cmd = """syst_fake=tf([1],[1 2 3]);
    syst_fake_dis=c2d(syst_fake,0.01);
    [output,t]=step(syst_fake_dis);
    out_len = length(output);
    input=1:out_len;
    input(:)=1;
    [num,den]=stmcb(output,input,0,2);
    sys_model=tf(num,den,0.01);"""
    out = matlab.eval(cmd)
    # Get matlab inputs
    #matlab_syst_fake = matlab.get('syst_fake')
    #matlab_syst_fake_dis = matlab.get('syst_fake_dis')
    matlab_output = matlab.get('output')
    matlab_t = matlab.get('t')
    matlab_out_len = matlab.get('out_len')
    matlab_input = matlab.get('input')
    # Get matlab output of stmcb
    matlab_num = matlab.get('num')
    matlab_den = matlab.get('den')
    #matlab_sys_model = matlab.get('sys_model')
     
    # Use same inputs on python stmcb
    [py_num, py_den] = stmcb(matlab_output,matlab_input,0,2)
     
    # Check if we get the same result
    print "\nMatlab num vs Python num"
    if np.allclose(matlab_num, py_num): # see if they are equal with relative tolerance=1e-05 and abs tol=1e-08
        print OKGREEN + "Coincident results!" + ENDC
    else:
        print FAIL +"Different results..." + ENDC
    print matlab_num
    print py_num
     
    print "\nMatlab den vs Python den"
    if np.allclose(matlab_den, py_den):
        print OKGREEN + "Coincident results!" + ENDC
    else:
        print FAIL +"Different results..." + ENDC
    print matlab_den
    print py_den
    
    
    
    