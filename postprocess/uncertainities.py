#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 08:59:04 2025

@author: av284181
"""

import numpy as np

def phase_uncertainty(z, delta_z):
    """
    Compute the uncertainty in the phase (argument) of a complex number a + bi,
    given uncertainties c (in a) and d (in b).
    
    Returns the uncertainty in radians.
    """
    a = np.real(z)
    b = np.imag(z)
    c = np.real(delta_z)
    d = np.imag(delta_z)
    
    numerator = 180/np.pi * np.sqrt((b * c)**2 + (a * d)**2)
    denominator = a**2 + b**2
    
    if denominator.any() == 0:
        raise ValueError("Undefined phase for zero complex number (a = b = 0).")
    
    delta_theta = numerator / denominator
    return delta_theta


def amplitude_uncertainty(z, delta_z):
    """
    Compute the uncertainty in the amplitude (modulus) of a complex number a + bi,
    given uncertainties c (in a) and d (in b), using NumPy.

    Returns the uncertainty in the amplitude.
    """
    a = np.real(z)
    b = np.imag(z)
    c = np.real(delta_z)
    d = np.imag(delta_z)
    
    numerator = np.sqrt((a * c)**2 + (b * d)**2)
    denominator = np.sqrt(a**2 + b**2)
    
    if denominator.any() == 0:
        raise ValueError("Undefined amplitude for zero complex number (a = b = 0).")

    delta_r = numerator / denominator
    return delta_r