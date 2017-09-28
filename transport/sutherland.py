"""
Sutherland Transport model.

    The file contains routines to compute transport properties using
    Sutherland's Law.
"""

import numpy as np


def computeviscosity(T, Cmu, S_ref):
    """Compute viscosity using sutherlands law."""
    mu = Cmu*(T*np.sqrt(T))/(T + S_ref)
    return mu


def computethermalconductivity(Pr, cp, mu):
    """Compute thermal conductivity."""
    kf = Pr*cp*mu
    return kf
