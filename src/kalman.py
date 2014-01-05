#!/usr/bin/env python
""" Computes the discrete Kalman filter """
from __future__ import division, print_function

__author__ = 'Richard Lindsley'
import numpy as np

class Kalman():
    def __init__(self, A, C, Qw, Qv, x00, P00):
        """ Initializes the discrete Kalman filter

        A: A(n-1) in x(n) = A(n-1)x(n-1) + w(n). A matrix, generally.
        C: C(n) in y(n) = C(n)x(n) + v(n). A matrix, generally.
        Qw: Autocorrelation matrix of w(n)
        Qv: Autocorrelation matrix of v(n)
        x00: Initial condition for x_hat(0|0)
        P00: Initial condition for P(0|0)

        Note that these values are stored as NumPy matrices to simplify the
        code below (due to multiple matrix multiplies needed)

        """
        self._A = np.matrix(A)
        self._C = np.matrix(C)
        self._Qw = np.matrix(Qw)
        self._Qv = np.matrix(Qv)
        self._x00 = np.matrix(x00)
        self._P00 = P00

        self._n = 1 # time index
        self._x_n1_n1 = self._x00 # x_hat(n-1|n-1), initialized to x_hat(0|0)
        self._P_n1_n1 = self._P00 # P(n-1|n-1), initialized to P(0|0)

        # Computed later
        self._x_n_n1 = None # x_hat(n|n-1)
        self._x_n_n = None # x_hat(n|n)
        self._P_n_n1 = None # P(n|n-1)
        self._P_n_n = None # P(n|n)
        self._K = None # K(n), the Kalman gain

    def predict(self):
        """ Predict x(n|n-1) """
        self._x_n_n1 = self._A * self._x_n1_n1
        self._P_n_n1 = self._A * self._P_n1_n1 * self._A.H + self._Qw
        self._K = self._P_n_n1 * self._C.H * (self._C * self._P_n_n1 *
                self._C.H + self._Qv).getI()
        return self._x_n_n1

    def update(self, y):
        """ Update x(n|n) with y(n) """
        self._x_n_n = self._x_n_n1 + self._K * (y - self._C * self._x_n_n1)
        i = np.eye(self._P_n_n1.shape[0])
        self._P_n_n = (i - self._K * self._C) * self._P_n_n1
        return self._x_n_n

    def timestep(self):
        """ Increment n and set variables accordingly """
        self._n += 1
        self._x_n1_n1 = self._x_n_n
        self._P_n1_n1 = self._P_n_n

    def print_status(self):
        #print('n={} P(n|n-1)={:0.4} K(n)={:0.4} P(n|n)={:0.4}'.format(self._n,
        print('----------------')
        print(' n={}'.format(self._n))
        print(' x(n|n)=\n{}'.format(self._x_n_n))
        print(' P(n|n-1)=\n{}'.format(self._P_n_n1))
        print(' K(n)=\n{}'.format(self._K))
        print(' P(n|n)=\n{}'.format(self._P_n_n))
