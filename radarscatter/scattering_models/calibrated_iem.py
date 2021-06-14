import numpy as np
import radarscatter
from radarscatter.scattering_models.fung_1992 import fung_1992


def calibrated_iem(theta_deg, epsilon_real, epsilon_imag, freq_GHz, pq, s_cm, error_threshold=1e-8):
    """
    Implements the calibrated IEM as outlined in Choker et al., 2017, based
    on the calibration of Baghdadi et al. in L, C, and X band.

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

    error_threshold : float
        Error threshold for truncating infinite series in IEM, defaults to 1e-8

    Outputs:
    --------
    sigma_dB
        Normalized backscatter in dB
    """
    theta_rad = theta_deg * np.pi/180

    # Calibrated Correlation Lengths
    if 1 <= freq_GHz <= 2:
        l_cm_hh = 2.6590 * theta_rad ** -1.4493 + 3.0484 * s_cm[i] * theta_rad ** -0.8044
        l_cm_vv = 5.8735 * theta_rad ** -1.0814 + 1.3015 * s_cm[i] * theta_rad ** -1.4498
    # C-band
    elif 4 <= freq_GHz <= 8:
        l_cm_hh = 0.162 + 3.006 * (np.sin(1.23 * theta_rad)) ** -1.494 * s_cm[i]
        l_cm_vv = 1.281 + 0.134 * (np.sin(0.19 * theta_rad)) ** -1.590 * s_cm[i]
    # X-band
    elif 8 <= freq_GHz <= 12:
        l_cm_hh = 18.102 * np.exp(-1.8908 * theta_rad) * s_cm[i] ** (0.7644 * np.exp(+0.2005 * theta_rad))
        l_cm_vv = 18.075 * np.exp(-2.1715 * theta_rad) * s_cm[i] ** (1.2594 * np.exp(-0.8308 * theta_rad))

    if pq == 'hh':
        l_cm = l_cm_hh
    elif pq == 'vv':
        l_cm = l_cm_vv
    elif pq == 'hv':
        raise ValueError("Cross-polarization not supported by 'calibrated_iem'")

    sigma_dB = fung_1992(theta_deg=theta_deg,
                         epsilon_real=epsilon_real,
                         epsilon_imag=epsilon_imag,
                         freq_GHz=freq_GHz,
                         pq=pq,
                         s_cm=s_cm,
                         l_cm=l_cm,
                         alpha=2,
                         use_transition_function=False,
                         error_threshold=error_threshold)
    return sigma_dB
