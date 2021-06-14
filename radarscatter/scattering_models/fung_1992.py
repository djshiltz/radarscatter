import numpy as np
import radarscatter
from radarscatter.miscellaneous.transition_function import fresnel_transition
from radarscatter.miscellaneous.trig import sind, cosd
from radarscatter.miscellaneous.power_spectrum import w_n_isotropic
import math


def fung_1992(theta_deg, epsilon_real, epsilon_imag, freq_GHz, pq, s_cm, l_cm, alpha, use_transition_function,
              error_threshold=1e-8):
    """
    Implements the original integral equation model published by Fung et al., 1992.

    The simplified single scattering equations are described in Ch. 3 of
    "Microwave Scattering and Emission Models for Users" by Fung et al. 2009.

    Inputs:
    -------
    theta_deg : float
        Incidence angle in degrees

    epsilon_real, epsilon_imag : float
        Real and (negative) imaginary parts of the dielectric ratio (soil to air)

    freq_GHz : float
        Radar frequency in GHz

    pq : str
        Polarization - either 'hh' or 'vv'

    s_cm : float
        RMS height of soil surface in cm

    l_cm : float
        Correlation length of soil surface in cm

    alpha : float
        Autocorrelation function shape parameter
        1 - exponential
        2 - Gaussian

    use_transition_function : bool
        Whether to use the transition function outlined in
        Chen et al., 2001

    error_threshold : float
        Error threshold for truncating infinite series, defaults to 1e-8

    Outputs:
    --------
    sigma_dB
        Normalized backscatter in dB
    """

    # Constants
    ep_r = epsilon_real - 1j * epsilon_imag
    mu_r = 1
    c = 2.998e+10  # speed of light in cm/s
    lambda_cm = c / (freq_GHz * 1e+9)
    k_radpcm = 2 * np.pi / lambda_cm
    sin = sind(theta_deg)
    cos = cosd(theta_deg)
    ks = k_radpcm * s_cm

    # Fresnel Reflection Coefficients
    R_h_theta, R_h_transition, R_v_theta, R_v_transition = fresnel_transition(epsilon=ep_r,
                                                                              theta_deg=theta_deg,
                                                                              k_radpcm=k_radpcm,
                                                                              s_cm=s_cm,
                                                                              l_cm=l_cm,
                                                                              alpha=alpha)
    if use_transition_function == True:
        R_h = R_h_transition
        R_v = R_v_transition
    elif use_transition_function == False:
        R_h = R_h_theta
        R_v = R_v_theta

    # Equation 3.1
    f_vv =  2 * R_v / cosd(theta_deg)
    f_hh = -2 * R_h / cosd(theta_deg)
    T_v =  1 + R_v
    T_vm = 1 - R_v
    T_h =  1 + R_h
    T_hm = 1 - R_h
    sq = np.sqrt(ep_r - sind(theta_deg)**2)
    F_vv =   (sin**2 / cos - sq / ep_r) * T_v**2 \
                - 2 * sin**2 * (1 / cos + 1 / sq) * T_v * T_vm \
                + (sin**2 / cos + ep_r * (1 + sin**2) / sq) * T_vm**2
    F_hh = -((sin**2 / cos - sq / mu_r) * T_h**2 \
                - 2 * sin**2 * (1 / cos + 1 / sq) * T_h * T_hm \
                + (sin**2 / cos + mu_r * (1 + sin**2) / sq) * T_hm**2)

    if pq == 'hh':
        f_pp = f_hh
        F_pp = F_hh
    elif pq == 'vv':
        f_pp = f_vv
        F_pp = F_vv
    elif pq == 'hv':
        raise ValueError("Cross-polarization not supported by 'fung_1992'")

    # Compute summation in Eq. 3.1
    summation = 0
    n = 1
    converged = False
    while not converged:
        I_ppn = (2 * ks * cos)**n * f_pp * np.exp(-(ks * cos)**2) + (ks * cos)**n * F_pp
        w_n = w_n_isotropic(k_r=2 * k_radpcm * sin,
                            n=n,
                            l_cm=l_cm,
                            alpha=alpha)
        term = np.abs(I_ppn)**2 * w_n / math.factorial(n)
        summation += term
        n += 1
        error = np.abs(term/summation)
        if np.all(error < error_threshold):
            converged = True
    sigma = k_radpcm**2 / (4 * np.pi) * np.exp(-2 * (ks * cos)**2) * summation
    sigma_dB = 10 * np.log10(sigma)
    return sigma_dB
