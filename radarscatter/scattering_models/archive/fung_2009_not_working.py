import numpy as np
import math
import radarscatter
from radarscatter.miscellaneous.trig import sind, cosd
from radarscatter.miscellaneous.transition_function import fresnel_transition
from radarscatter.miscellaneous.power_spectrum import w_n_isotropic


validation = True

def fung_2009(freq_GHz, pq, epsilon_real, epsilon_imag, theta_deg, s_cm, l_cm, alpha):
    """
    This is the IEM-B model described in Ch. 4 of "Microwave Scattering and Emission
    Models for Users" by Fung and Chen.  It includes the transition function described
    on page 164.

    This code assumes mu_r = 1 (non-magnetic soil)

    Inputs:
    -------
    freq_GHz : float
        radar frequency in GHz

    pq : str
        radar polarization, either 'hh' or 'vv'

    epsilon_real, epsilon_imag : float
        real and imaginary parts of the (relative) dielectric constant of the soil

    theta_deg : float
        incidence angle in degrees

    s_cm, l_cm : float
        RMS height and correlation length of the surface in cm

    alpha : float
        shape parameter of the autocorrelation function, 1: exponential, 2: gaussian

    Outputs:
    --------
    sigma_dB : float
        normalized radar cross section in dB
    """

    c_cmps = 2.998e+10  #cm/s
    freq_Hz = freq_GHz * 1e-9
    lambda_cm = c_cmps / freq_Hz
    k_radpcm = 2 * np.pi / lambda_cm
    kz = k_radpcm * cosd(theta_deg)
    epsilon = epsilon_real - 1j * epsilon_imag
    sq = np.sqrt(epsilon - sind(theta_deg)**2)

    R_h, R_ht, R_v, R_vt = fresnel_transition(epsilon=epsilon,
                                              theta_deg=theta_deg,
                                              k_radpcm=k_radpcm,
                                              s_cm=s_cm,
                                              l_cm=l_cm,
                                              alpha=alpha)
    f_vv = 2 * R_vt / cosd(theta_deg)
    f_hh = -2 * R_ht / cosd(theta_deg)
    F_vv1 =  4 * k_radpcm / sq * ((1 - R_vt)**2 * epsilon * cosd(theta_deg) \
                + (1 - R_vt) * (1 + R_vt) * sind(theta_deg)**2 * (sq - cosd(theta_deg)) \
                - (1 + R_vt)**2 * (cosd(theta_deg) + sind(theta_deg)**2 / (2 * epsilon) * (sq - cosd(theta_deg))))
    F_hh1 = -4 * k_radpcm / sq * ((1 - R_ht)**2 * cosd(theta_deg) \
                + (1 - R_ht) * (1 + R_ht) * sind(theta_deg)**2 * (sq - cosd(theta_deg)) \
                - (1 + R_ht)**2 * (epsilon * cosd(theta_deg) + sind(theta_deg)**2 / 2 * (sq - cosd(theta_deg))))
    F_vv2 =  4 * k_radpcm * sind(theta_deg)**2 * ((1 - R_vt)**2 * (1 + epsilon * cosd(theta_deg) / sq) \
                - (1 - R_vt) * (1 + R_vt) * (3 + cosd(theta_deg) / sq) \
                + (1 + R_vt)**2 * (1 + 1 / (2 * epsilon) + epsilon * cosd(theta_deg) / (2 * sq)))
    F_hh2 = -4 * k_radpcm * sind(theta_deg)**2 * ((1 - R_ht)**2 * (1 + cosd(theta_deg) / sq) \
                - (1 - R_ht) * (1 + R_ht) * (3 + cosd(theta_deg) / sq) \
                + (1 + R_ht)**2 * (1 + 1/2 + cosd(theta_deg) / (2 * sq)))

    f_pp =  [f_hh, f_vv]
    F_pp1 = [F_hh1, F_vv1]
    F_pp2 = [F_hh2, F_vv2]

    if pq == 'hh':
        f_pp = f_hh
        F_pp1 = F_hh1
        F_pp2 = F_hh2
    elif pq == 'vv':
        f_pp = f_vv
        F_pp1 = F_vv1
        F_pp2 = F_vv2
    elif pq == 'hv':
        raise ValueError("Cross-polarization not supported by 'fung_2009'")

    threshold = 1e-8
    sum = np.float64(0)
    n = 2
    converged = False
    while not converged:
        term = np.abs((2 * kz * s_cm)**n * f_pp + s_cm/4 * F_pp1 * (2 * kz * s_cm)**(n - 1))**2 \
                   * w_n_isotropic(k_r=2 * k_radpcm * sind(theta_deg),
                                   n=n,
                                   l_cm=l_cm,
                                   alpha=alpha) / math.factorial(n)
        sum += term.astype(float)
        n += 1
        error = np.abs(term / sum)
        if error < threshold:
            converged = True

    sigma_0 = k_radpcm**2 / (4 * np.pi) * np.exp(-4 * kz**2 * s_cm**2) \
                  * (np.abs((2 * kz * s_cm) * f_pp + s_cm/4 * (F_pp1 + F_pp2))**2 \
                  * w_n_isotropic(k_r=2 * k_radpcm * sind(theta_deg),
                                  n=1,
                                  l_cm=l_cm,
                                  alpha=alpha) \
                  + sum)
    sigma_dB = 10 * np.log10(sigma_0)
    return sigma_dB




if validation:
    import matplotlib.pyplot as plt
    plt.style.use('../../custom.mplstyle')

    # Figures 4.2 and 4.3
    theta_deg = np.linspace(start=0.1, stop=70, num=1000)
    sigma_dB = np.zeros_like(theta_deg)
    freq_GHz = 5
    l_cm = 5
    alpha = 1
    epsilon_real = [6, 36]
    epsilon_imag = [0.2, 1.2]
    s_cm = [0.2, 0.3, 0.4, 0.6, 0.8]
    colors = ['red', 'olive', 'cyan', 'brown', 'magenta']

    fig1 = plt.figure(num=1, figsize=(15, 15))
    ax1 = fig1.add_subplot(221)
    for j, s_cm_j in enumerate(s_cm):

        for i, theta_deg_i in enumerate(theta_deg):
            sigma_dB[i] = fung_2009(freq_GHz=freq_GHz,
                                    pq='vv',
                                    epsilon_real=epsilon_real[0],
                                    epsilon_imag=epsilon_imag[0],
                                    theta_deg=theta_deg_i,
                                    s_cm=s_cm_j,
                                    l_cm=l_cm,
                                    alpha=alpha)
        ax1.plot(theta_deg, sigma_dB, color=colors[j])
    ax1.set_xlabel(r'$\theta$ [$^\circ$]')



