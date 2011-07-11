"""Python code to calculate the abscissas and weights needed for
Legendre-Gauss quadrature.  These routines are not optimized but work fast
enough for l ~ 500 that I am not worried about speeding them up."""

# $Id: legendre_gauss.py,v 1.2 2011-07-08 05:54:57 copi Exp $

import scipy.special as sf
import scipy.optimize as opt
import numpy as np

def Pl_func_theta (theta, l) :
    """Compute Pl(cos theta).  This routine is NOT optimized."""
    return Pl_func (np.cos(theta), l)

def Pl_func (x, l) :
    """Compute Pl(x).  This routine is NOT optimized."""
    (pl, dpl) = sf.lpn (l, x)
    return pl[l]

def abscissas (l, tol=1.0e-12) :
    """Return the abscissas, the positive zeros of P_l(cos theta).  The
    values of cos(theta) are returned.  Since the P_l are (anti)symmetric
    for l (odd) even all the zeros can easily be determined from this list.
    For l odd cos(theta)=0 is a zero and is included in the list."""
    Lmax = int((l + 1) / 2)
    Plzero = np.zeros(Lmax)
    for j in range (Lmax) :
        tmin = j/(l+0.5)*np.pi
        tmax = (j+1)/(l+0.5)*np.pi
        Plzero[Lmax-1-j] = opt.zeros.bisect (Pl_func_theta, tmin, tmax, args=(l),
                                             xtol=tol)
    Plzero = np.cos(Plzero)
    if l%2 == 1 : Plzero[0] = 0
    return Plzero

def weights (l, Plzeros) :
    """Return the weights for Legendre-Gauss quadrature.  The zeros should
    be calculated using abscissas().  For even l: 2*weights.sum() == 2; for
    odd l: 2*weights.sum() - weights[0] == 2; since +/- of the zeros
    provided are zeros."""
    w = np.zeros(Plzeros.size)
    for j in range(w.size) :
        w[j] = 2 * (1-Plzeros[j]**2) / ((l+1)*Pl_func(Plzeros[j], l+1))**2
    return w

def abscissas_weights (l, tol=1.0e-12) :
    """Return a tuple containing the abscissas (zeros of the Legendre
    polynomial of order l) and the weights for Legendre-Gauss quadrature.
    See abscissas() and weights() for more details.
    
    This routine returns the FULL list of zeros and weights, that is, the
    zeros at positive and negative arguments to the Legendre polynomial
    unlike abscissas() and weights().
    """
    z = abscissas (l, tol)
    w = weights (l, z)
    if l%2 == 0 : # Even
    	Z = np.r_ [ -z[::-1], z ]
    	W = np.r_ [  w[::-1], w ]
    else : # Odd, do not double count the "0"
    	Z = np.r_ [ -z[:0:-1], z ]
    	W = np.r_ [  w[:0:-1], w ]
    return (Z, W)
