#!/usr/bin/env python
""" LMS adaptive filters """
from __future__ import print_function, division
__author__ = 'Richard Lindsley'

import sys
import numpy as np
import matplotlib.pyplot as plt
import convm

def lms(x, d, p, mu, w0=None):
    """ Finds the LMS adaptive filter w based on the inputs

    Adapted from the Hayes MATLAB code.

    Inputs:
        x: input data to the filter, Nx1 vector
        d: desired output signal, Nx1 vector
        p: filter order
        mu: step size, scalar or Nx1 vector
        w0: initial guess (optional)

    Outputs:
        w: filter coefficient matrix: w[n,p] where 
            n is the time index and 
            p is the filter tap index
        e: error sequence e(n): Nx1 vector

    """
    # Initialize
    N = len(x)
    X = convm.convm(x, p)
    if not w0:
        w0 = np.zeros(p)
    # Promote mu to be a vector if it's not already
    if not isinstance(mu, np.ndarray):
        mu = np.ones(N) * mu
    w = np.zeros((N,p))
    e = np.zeros(N)
    w[0,:] = w0

    # Run the filter
    for n in xrange(N-1):
        y = np.dot(w[n,:], X[n,:])
        e[n] = d[n] - y
        w[n+1,:] = w[n,:] + mu[n] * e[n] * X[n,:].conj()

    return w, e

def nlms(x, d, p, beta, w0=None, verbose=False):
    """ Finds the NLMS adaptive filter w based on the inputs

    Adapted from the Hayes MATLAB code.

    Inputs:
        x: input data to the filter, Nx1 vector
        d: desired output signal, Nx1 vector
        p: filter order
        beta: normalized step size (0 < beta < 2)
        w0: initial guess (optional)
        verbose: verbosity flag

    Outputs:
        w: filter coefficient matrix: w[n,p] where 
            n is the time index and 
            p is the filter tap index
        e: error sequence e(n): Nx1 vector

    """
    # Initialize
    N = len(x)
    eps = 1e-4 # epsilon for denominator
    X = convm.convm(x, p)
    if not w0:
        w0 = np.zeros(p)
    w = np.zeros((N,p))
    e = np.zeros(N)
    w[0,:] = w0
    if verbose:
        print(' ')

    # Run the filter
    for n in xrange(N-1):
        if verbose and n % 10000 == 0:
            print(' {}/{}\r'.format(n,N-1), end='')
            sys.stdout.flush()
        y = np.dot(w[n,:], X[n,:])
        e[n] = d[n] - y
        norm = eps + np.dot(X[n,:], X[n,:].T.conj())
        w[n+1,:] = w[n,:] + beta / norm * e[n] * X[n,:].conj()

    if verbose:
        print(' ')
    return w, e

def plot_traj(w, w_true, fname=None):
    """ Plot the trajectory of the filter coefficients

    Inputs:
        w: matrix of filter weights versus time
        w_true: vector of true filter weights 
        fname: what filename to export the figure to. If None, then doesn't
        export

    """
    plt.figure()
    n = np.arange(w.shape[0])
    n_ones = np.ones(w.shape[0])
    plt.hold(True)
    # NOTE: This construction places a limit of 4 on the filter order
    plt_colors = ['b', 'r', 'g', 'k']
    for p in xrange(w.shape[1]):
        plt.plot(n, w[:,p], '{}-'.format(plt_colors[p]), label='w({})'.format(p))
        plt.plot(n, w_true[p] * n_ones, '{}--'.format(plt_colors[p]))
    plt.xlabel('Iteration')
    plt.ylabel('Coefficients')
    plt.legend()
    if fname:
        plt.savefig(fname)
        plt.close()

def plot_error(e, fname=None):
    """ Plot the squared error versus time

    Inputs:
        e: vector of the error versus time
        fname: what filename to export the figure to. If None, then doesn't
        export

    """
    plt.figure()
    n = np.arange(len(e))
    e2 = np.power(e, 2)
    plt.plot(n, e2)
    #plt.semilogy(n, e2)
    plt.xlabel('Iteration')
    plt.ylabel('Squared error')
    if fname:
        plt.savefig(fname)
        plt.close()
