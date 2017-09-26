"""
Calorically perfect gas model.

    The file contains routines to compute thermodynamic properties using the
    calorically perfect gas model.
"""


def computedensity(P, T, R):
    """Compute density using ideal gas law."""
    rho = P/(R*T)
    return rho


def computeinternalenergy(P, gamma, rho):
    """Compute interal engery using cpg (constant gamma)."""
    ein = P/((gamma-1.0)*rho)
    return ein


def computetotalenergy(ein, u, v):
    """Compute total engergy from internal and kinetic engery."""
    Et = ein + 0.5*(u**2.0 + v**2.0)
    return Et


def computecp(gamma, R):
    """Compute specific hear for constant pressure for cpg (constant gamma)."""
    cp = (gamma*R)/(gamma-1.0)
    return cp


def computecv(gamma, R):
    """Compute specific hear for constant volume for cpg (constant gamma)."""
    cv = R/(gamma-1.0)
    return cv
