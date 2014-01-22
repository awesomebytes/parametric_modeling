
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

syst_fake = tf(np.array(np.hstack((1.))), np.array(np.hstack((1., 2., 3.))))
syst_fake_dis = c2d(syst_fake, 0.01)
[output, t] = plt.step(syst_fake_dis)
plt.plot(output)
out_len = length(output)
input = np.arange(1., (out_len)+1)
input[:] = 1.
[num, den] = stmcb_test(output, input, 0., 2.)
#% sys_model=tf(num,den,0.01)
#% step(sys_model)
#% hold on
#% step(syst_fake)