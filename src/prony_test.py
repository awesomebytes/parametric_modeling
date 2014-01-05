#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 5 17:11:46 2014

@author: Sammy Pfeiffer
This file is a test for using prony python implementation versus MATLAB one.
"""

import numpy as np
from prony import prony
from scipy.signal import butter, impulse
from impz import impz

# %   % Example:
# %   %   Fit an IIR model to an impulse response of a lowpass filter.
# %
# %   [b,a] = butter(4,0.2);
b, a = butter(4, 0.2)
# MATLAB output
# Python output

# b = 0.0048    0.0193    0.0289    0.0193    0.0048
#[ 0.00482434  0.01929737  0.02894606  0.01929737  0.00482434]

# a = 1.0000   -2.3695    2.3140   -1.0547    0.1874
#[ 1.         -2.36951301  2.31398841 -1.05466541  0.18737949]
print b
print a

# %   impulseResp = impz(b,a);                % obtain impulse response
# 0.00482434335771622    0.0307287177680857    0.0905946819548826    0.167944821844737    0.224641271344028    0.233457187867600    0.193512552162805    0.123765243571016    0.0496036031380585    -0.00850905187491667    -0.0406738350178057    -0.0475631979469677    -0.0368517338223919    -0.0185628385243508    -0.00125221912683403    0.0100331628527336    0.0139990059845164    0.0121118272327115    0.00712186446378598    0.00173298095479121    -0.00222279239538353    -0.00403535730812247    -0.00392509206552861    -0.00263181398008669    -0.000992945935699223    0.000353673136269240    0.00109549708726230    0.00122332129708734    0.000922772684652072    0.000444882357824948    3.79019631817530e-06    -0.000276480597968442    -0.000367601488225979    -0.000310628048735716    -0.000177716344133915    -3.82012617315314e-05    6.19874979750476e-05    0.000106051505667546    0.000100862919097785    6.61282460408185e-05    2.35297812645709e-05    -1.07611182300523e-05    -2.91027199865281e-05
impulseResp = impz(b,a)
print impulseResp
# I used a couple of examples... matlab gave 43 as length we give 64
# python [ 0.00482434  0.03072872  0.09059468  0.16794482  0.22464127]
# matlab  0.00482434335771622    0.0307287177680857    0.0905946819548826    0.167944821844737    0.224641271344028

# %   denOrder=4; numOrder=4;                 % system function of order 4
denOrder = 4
numOrder = 4

# %   [Num,Den]=prony(impulseResp,numOrder,denOrder);
[Num,Den, err]=prony(impulseResp,numOrder,denOrder)
# MATLAB
# Num =0.0048    0.0193    0.0289    0.0193    0.0048
# Den =1.0000   -2.3695    2.3140   -1.0547    0.1874
# Python
# [ 0.00482434  0.01929737  0.02894606  0.01929737  0.00482434]
# [ 1.         -2.36951301  2.31398841 -1.05466541  0.18737949]

print Num
print Den

# %   subplot(211);                           % impulse response and input
# %   stem(impz(Num,Den,length(impulseResp)));   
# %   title('Impulse Response with Prony Design');
# %   subplot(212);
# %   stem(impulseResp); title('Input Impulse Response');


