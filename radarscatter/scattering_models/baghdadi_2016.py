import numpy as np

def baghdadi_2016(ks : float,
                  pq : str,
                  theta_deg : float,
                  mv_pct : float):
    """
    Computes radar backscatter based on the empirical model from Baghdadi et al., 2016
    """

    if pq == 'hh':
        # eq. (5)
        sigma    = 10 ** -1.287 * cos(theta_deg) ** 1.227 \
                   * 10 ** (0.009 * cot(theta_deg) * mv_pct) \
                   * ks ** (0.86 * sin(theta_deg))

    elif pq == 'vv':
        # eq. (6)
        sigma    = 10 ** -1.138 * cos(theta_deg) ** 1.528 \
                   * 10 ** (0.008 * cot(theta_deg) * mv_pct) \
                   * ks ** (0.71 * sin(theta_deg))

    elif pq == 'hv':
        # eq. (7)
        sigma    = 10 ** -2.325 * cos(theta_deg) ** -0.01 \
                   * 10 ** (0.011 * cot(theta_deg) * mv_pct) \
                   * ks ** (0.44 * sin(theta_deg))

    sigma_dB = 10 * np.log10(sigma)

    return sigma_dB


def cos(x_deg):
    return np.cos(x_deg * np.pi / 180)


def cot(x_deg):
    return 1/np.tan(x_deg * np.pi / 180)


def sin(x_deg):
    return np.sin(x_deg * np.pi / 180)
