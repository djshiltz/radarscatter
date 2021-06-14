import numpy as np


def dubois_1995(ks : float,
                lambda_cm : float,
                pq : str,
                theta_deg : float,
                epsilon_real : float,):
    """
    Computes radar backscatter based on the empirical model from Dubois et al., 1995

    Note: there were multiple typos in the original publication.  The equations given
    here are the correct ones according to Ulaby and Long, 2015 and Baghdadi et al., 2016

    Ulaby and Long, 2015 refer to this model as the "SMART" model
        - Soil Moisture Assessment Radar Technique
    """

    if pq == 'hh':
        # eq. (1a)
        sigma = 10 ** -2.75 * cos(theta_deg) ** 1.5 / sin(theta_deg) ** 5 \
              * 10 ** (0.028 * epsilon_real * tan(theta_deg)) \
              * (ks * sin(theta_deg)) ** 1.4 \
              * lambda_cm ** 0.7

    elif pq == 'vv':
        # eq. (1b)
        sigma = 10 ** -2.35 * cos(theta_deg) ** 3 / sin(theta_deg) ** 3 \
              * 10 ** (0.046 * epsilon_real * tan(theta_deg)) \
              * (ks * sin(theta_deg)) ** 1.1 \
              * lambda_cm ** 0.7

    elif pq == 'hv':
        raise ValueError("Cross-polarization not supported by 'dubois_1995'")

    sigma_dB = 10 * np.log10(sigma)

    return sigma_dB


def cos(x_deg):
    return np.cos(x_deg * np.pi / 180)


def sin(x_deg):
    return np.sin(x_deg * np.pi / 180)


def tan(x_deg):
    return np.tan(x_deg * np.pi/180)