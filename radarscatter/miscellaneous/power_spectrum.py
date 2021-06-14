import numpy as np


def w_n_isotropic(k_r, n, l_cm, alpha):
    """
    Computes the Hankel transform of the autocorrelation function (assumed azimuthally
    symmetric) to the nth power

    These formulas are found in the paper: "Sensitivity of Bistatic
    Scattering to Soil Moisture and Surface Roughness of Bare Soils" by
    Brogioni et al. 2010.

    In this paper, the 2D Fourier transform is normalized by the factor 2pi, so
    this function returns 2pi times the expressions listed in the paper so
    that the expression is consistent with other references on surface
    scattering.

    Inputs:
    -------
    k_r : float
        wavenumber (radians per cm) in the radial direction

    n : int
        power to which ACF is raised

    l_cm : float
        correlation length in cm

    alpha : int
        shape parameter of autocorrelation function: 1=exponential, 2=gaussian
    """
    if alpha not in [1, 2]:
        raise ValueError("Shape parameter 'alpha' must be in [1, 2]")

    if alpha == 1:
        return 2 * np.pi * (l_cm/n)**2 * (1 + k_r**2 * l_cm**2 / n**2)**(-1.5)

    if alpha == 2:
        return 2 * np.pi * l_cm**2 / (2*n) * np.exp(-k_r**2 * l_cm**2 / (4*n))
