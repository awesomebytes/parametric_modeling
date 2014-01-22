#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Jan 22 22:26 2014

@author: Sammy Pfeiffer

Test our python arcov vs the MatlabTM armcov

MATLAB testing code from:
http://www.music.mcgill.ca/~gary/307/matlab/predict.m
"""
from matlabpipe import MatlabPipe
import numpy as np
from armcov import armcov

OKGREEN = '\033[92m'
FAIL = '\033[91m'
ENDC = '\033[0m'

if __name__ == '__main__':
    matlab = MatlabPipe(matlab_version='2013a')
    matlab.open()
    
    cmd = """% Signal parameters.
        fs = 11025;
        T = 1/fs;
        N = 40000;              % signal length
        t60 = [1.0 0.8 0.6 0.7 0.7];    % T60 decay constants
        tau = -t60/log(0.001);  % 1/e decay constants
        f = [220 660 1100 4000 6000];     % sinusoid frequencies (Hz)
        phi = [pi/4 pi/2 pi 0 0];   % sinusoid phase offsets
        p = length(f);          % number of resonances
        p = 3;
        
        t = [0:N]*T;            % time vector
        x = 0;                  % initialize our input signal
        for n = 1:p,
          x = x + exp(-t/tau(n)).*sin(2*pi*f(n)*t+phi(n));
        end
        x = 0.99*x/max(abs(x));
        
        % Do linear prediction using covariance method.
        [a e] = armcov(x, 2*p);"""
    out = matlab.eval(cmd)
    # Get matlab inputs
    matlab_x = matlab.get('x')
    matlab_p = matlab.get('p')
    # Get matlab output of armcov
    matlab_a = matlab.get('a')
    matlab_e = matlab.get('e')

     
    # Use same inputs on python armcov
    [py_a, py_e] = armcov(matlab_x, 2*matlab_p)
     
    # Check if we get the same result
    print "\nMatlab a vs Python a"
    if np.allclose(matlab_a, py_a): # see if they are equal with relative tolerance=1e-05 and abs tol=1e-08
        print OKGREEN + "Coincident results!" + ENDC
    else:
        print FAIL +"Different results..." + ENDC
    print matlab_a
    print py_a
     
    
    print "\nMatlab e vs Python e"
    if np.allclose(matlab_e, py_e): # see if they are equal with relative tolerance=1e-05 and abs tol=1e-08
        print OKGREEN + "Coincident results!" + ENDC
    else:
        print FAIL +"Different results..." + ENDC
    print matlab_e
    print py_e
    