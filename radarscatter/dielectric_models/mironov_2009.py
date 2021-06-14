import numpy as np


def mironov_2009(clay_pct, mv_pct, freq_GHz):
    """
    Computes the dielectric constant (DC) and loss factor (LF)
    using the mineralogy-based soil dielectric model (MBSDM) described in
    the paper "Physically and Mineralogically Based Spectroscopic
    Dielectric Model for Moist Soils," by Mironov et. al 2009.

    Inputs:
    -------
    clay_pct : float
        clay mass fraction (0-100)[%]

    mv_pct : float
        volumetric moisture content (0-100)[%]

    freq_GHz : float
         frequency (GHz)

    Outputs:
    --------
    ep_r : complex
        complex (relative) dielectric constant of soil/water mixture (unitless)
    """

    freq_Hz = freq_GHz * 1e+9
    mv_dec = mv_pct/100

    # Empirical Curves (equations 17-25)
    sigma = {'b': 0.3112 + 0.467e-2 * clay_pct, \
             'u': 0.3631 + 1.217e-2 * clay_pct}

    e_0 = {'b': 79.8 - 85.4e-2 * clay_pct + 32.7e-4 * clay_pct ** 2, \
           'u': 100}  # bottom left of p. 2064

    tau = {'b': 1.062e-11 + 3.45e-12 * 1e-2 * clay_pct, \
           'u': 8.5e-12}  # bottom left of p. 2064

    e_inf = 4.9  # see paragraph below eq. 16
    e_0f = 8.854e-12  # F/m , top left of p. 2061

    # Dielectric properties of bound and unbound (free) water

    def DC_water_fcn(e_0, tau):
        # Equation 16a
        return e_inf + (e_0 - e_inf) / (1 + (2 * np.pi * freq_Hz * tau) ** 2)

    def LF_water_fcn(e_0, tau, sigma):
        # Equation 16b
        return (e_0 - e_inf) / (1 + (2 * np.pi * freq_Hz * tau) ** 2) * 2 * np.pi * freq_Hz * tau \
               + sigma / (2 * np.pi * e_0f * freq_Hz)

    DC_water = {'b': DC_water_fcn(e_0['b'], tau['b']), \
                'u': DC_water_fcn(e_0['u'], tau['u'])}

    LF_water = {'b': LF_water_fcn(e_0['b'], tau['b'], sigma['b']), \
                'u': LF_water_fcn(e_0['u'], tau['u'], sigma['u'])}

    # Refractive index (RI) and normalized attenuation coefficient (NAC)

    # eq. 14
    def n_fcn(DC, LF):
        return np.sqrt(np.sqrt(DC ** 2 + LF ** 2) + DC) / np.sqrt(2)

    # eq. 15
    def k_fcn(DC, LF):
        return np.sqrt(np.sqrt(DC ** 2 + LF ** 2) - DC) / np.sqrt(2)

    # eq. 17
    n = {'d': 1.634 - 0.539e-2 * clay_pct + 0.2748e-4 * clay_pct ** 2, \
         'b': n_fcn(DC_water['b'], LF_water['b']), \
         'u': n_fcn(DC_water['u'], LF_water['u'])}

    # eq. 18
    k = {'d': 0.03952 - 0.04038e-2 * clay_pct, \
         'b': k_fcn(DC_water['b'], LF_water['b']), \
         'u': k_fcn(DC_water['u'], LF_water['u'])}

    # Maximum bound water fraction (MBWF)

    # eq. 18
    m_vt = 0.02863 + 0.30673e-2 * clay_pct

    LEQ = (mv_dec <= m_vt)
    LEQ = int(LEQ)

    # eq. 12
    n_m = LEQ * (n['d'] + (n['b'] - 1) * mv_dec) \
          + (1 - LEQ) * (n['d'] + (n['b'] - 1) * m_vt + (n['u'] - 1) * (mv_dec - m_vt))

    # eq. 13
    k_m = LEQ * (k['d'] + k['b'] * mv_dec) \
          + (1 - LEQ) * (k['d'] + k['b'] * m_vt + k['u'] * (mv_dec - m_vt))

    # Return Result
    # eq. 11
    DC_m = n_m ** 2 - k_m ** 2
    LF_m = 2 * n_m * k_m

    return DC_m - 1j * LF_m
