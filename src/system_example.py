#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 12:30:46 2013

@author: Sammy Pfeiffer
"""
from numpy import allclose             # Grab all of the NumPy functions
import numpy as np
from matplotlib.pyplot import plot, show # Grab MATLAB plotting functions
from control.matlab import tf, c2d, step    # MATLAB-like functions
#from scipy.signal import cont2discrete as c2d
#from scipy.signal import step
from stmbc import stmcb


# MATLAB:
# syst_fake=tf([1],[1 2 3]) # from matlab to python, commas must be used to separate elements
# >> syst_fake
# syst_fake =
#  
#         1
#   -------------
#   s^2 + 2 s + 3
#  
# Continuous-time transfer function.

syst_fake = tf([1],[1, 2, 3])
print syst_fake
# Python:
#       1
# -------------
# s^2 + 2 s + 3


# MATLAB:
# zero order hold by default
# syst_fake_dis=c2d(syst_fake,0.01)
# syst_fake_dis =
#  
#   4.967e-05 z + 4.934e-05
#   -----------------------
#    z^2 - 1.98 z + 0.9802
#  
# Sample time: 0.01 seconds
# Discrete-time transfer function.

#Test if we generated a system with the same behaviour:
# MATLAB:
# [output, t] = step(syst_fake)
# plot(t, output)
# Python
# output, t = step([syst_fake.num[0][0][0], syst_fake.den[0][0]])
# plot(output, t) # MATLAB wants the inverse order to plot
# show()


syst_fake_dis = c2d(syst_fake, 0.01, method='zoh')
print syst_fake_dis
# OLD scipy.signal.cont2discrete way!
#syst_fake_dis = c2d([syst_fake.num[0][0][0], syst_fake.den[0][0]], 0.01)
#print_c2d_matlablike(syst_fake_dis)

# Python:
# zero order hold by default
# 4.96670866519e-05 z 4.93370716897e-05
# -----------------------
#  z^2 1.0 z^2 -1.97990166083 z 0.980198673307
# 
# Sample time: 0.01 seconds
# Discrete-time transfer function.

# MATLAB:
# [output,t]=step(syst_fake_dis)
### OLD scipy.signal.step way!
###output, t = step(syst_fake_dis[0][0], syst_fake_dis[1])# ,T=200) # MATLAB calculates the T size automatically based on black magic different than python one, see [MATLAB]/toolbox/shared/controllib/engine/@DynamicSystem/step.m:
output, t = step(syst_fake_dis)
# MATLAB:
# plot(output)
plot(output[0])
show()

# out_len = len(output)
out_len = len(output[0])

# input=1:650;
# input(:)=1;
input_ = np.ones(out_len)
# [num,den]=stmcb(output,input,0,2)
[num,den]=stmcb(output[0],input_,0,2)
# sys_model=tf(num,den,0.01)
# step(sys_model)
# hold on
# step(syst_fake)