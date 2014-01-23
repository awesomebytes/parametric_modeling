parametric_modeling
===================

Parametric modelling functions that can be found in MATLAB... now in Python!
(AR, ARMA, frequency response modeling)
The idea is to have all the functions listed at: http://www.mathworks.es/es/help/signal/parametric-modeling.html

Right now we have:

* **arburg**	Autoregressive (AR) all-pole model parameters estimated using Burg method
  * Working thanks to spectrum arburg (http://thomas-cokelaer.info/software/spectrum/html/user/ref_param.html#spectrum.burg.arburg)

* **arcov**	Estimate AR model parameters using covariance method
  * Working thanks to spectrum arcovar (http://thomas-cokelaer.info/software/spectrum/html/user/ref_psd_other.html#spectrum.covar.arcovar)

* **armcov**	Estimate AR model parameters using modified covariance method
  * Working thanks to spectrum modcovar (http://thomas-cokelaer.info/software/spectrum/html/user/ref_psd_other.html#spectrum.modcovar.modcovar)

* **aryule**	Estimate autoregressive (AR) all-pole model using Yule-Walker method
  * Working thanks to spectrum aryule (http://thomas-cokelaer.info/software/spectrum/html/user/ref_param.html#spectrum.yulewalker.aryule)

* **invfreqs**	Identify continuous-time filter parameters from frequency response data
  * Working thanks to a file found at http://projects.scipy.org/scipy/attachment/ticket/393/invfreq.py

* **invfreqz**	Identify discrete-time filter parameters from frequency response data
  * Not working :( whenever polystab function in spectrum is implemented this file can be finalized easily (it's nearly done!)

* **prony**	 Prony method for filter design
  * Working thanks to... me!

* **stmcb**	Compute linear model using Steiglitz-McBride iteration
  * Working thanks to... me!

Also I needed to implement:

* **convmtx**    Convolution matrix
  * Working thanks to... me!


===================

Developed using **Python 2.7, Eclipse + PyDev, MATLAB R2013a 64 bits**

* Using **numpy 1.6.1** (from ubuntu debs, using Ubuntu 12.04 64 bit)
```
sudo apt-get install numpy
```

* Using **scipy 0.13.1** (upgraded it from Ubuntu debs as I needed newer functions and some bugfixes)
```
sudo pip install scipy
```

* Using **spectrum 0.5.6** https://pypi.python.org/pypi/spectrum
```
sudo pip install spectrum 
```


I wanted to use/integrate in: **python_control**, using a checkout of the most updated branch, can be found at:
https://github.com/awesomebytes/python-control-code

But I don't have time for now.


* For testing purposes using **matlabpipe** from **python-mlabwrap**

This great library lets you **execute MATLAB code from Python** (easily! check the tests xxxxx_matlab_vs_python.py).

https://github.com/awesomebytes/python-mlabwrap

This is a fork from the work of Filipe Fernandes which he forked from http://code.google.com/p/danapeerlab/
```
git clone https://github.com/awesomebytes/python-mlabwrap
cd python-mlwrap
sudo python setup.py install
```

Thanks to this tool I can test the behaviour of my functions versus the MATLAB functions
and also give the same input data without needing to implement lots of other functions.

* Also another incredible tool is **libermate** http://libermate.sourceforge.net/

**A MATLAB to Python code translator** which doesn't do all the job but helps

I've rehosted it now in Github https://github.com/awesomebytes/libermate to give it more visibility and 
add it's README to it and my advice on how to use it if you ever want to translate a MATLAB script to Python.


