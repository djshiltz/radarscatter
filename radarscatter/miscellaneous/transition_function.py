import numpy as np
import math
from radarscatter.miscellaneous.fresnel_equations import fresnel
from radarscatter.miscellaneous.power_spectrum import w_n_isotropic
from radarscatter.miscellaneous.trig import sind, cosd


validation = False

def fresnel_transition(epsilon, theta_deg, k_radpcm, s_cm, l_cm, alpha):
    """
    Implements the transition function for Fresnel reflectance coefficients
    as described in "Microwave Scattering and Emission Models
    for Users" by Fung et al. 2009.

    These formulas are given on page 164

    I assume mu_r = 1 (non-magnetic soil)
    """

    R_h, R_v = fresnel(epsilon=epsilon,
                       theta_deg=theta_deg)
    R_h0, R_v0 = fresnel(epsilon=epsilon,
                         theta_deg=0)

    F_t = 8 * R_v0 ** 2 * sind(theta_deg) ** 2 \
          * (cosd(theta_deg) + np.sqrt(epsilon - sind(theta_deg) ** 2) \
          / (cosd(theta_deg) * np.sqrt(epsilon - sind(theta_deg) ** 2)))

    S_t0 = np.abs(1 + 8 * R_v0 / (F_t * cosd(theta_deg))) ** (-2)

    threshold = 1e-8

    sum_numerator = np.float64(0)
    n = 1
    converged = False
    while not converged:
        term = (k_radpcm * s_cm * cosd(theta_deg)) ** (2 * n) / math.factorial(n) \
               * w_n_isotropic(k_r=2 * k_radpcm * sind(theta_deg),
                               n=n,
                               l_cm=l_cm,
                               alpha=alpha)
        sum_numerator += term.astype(float)
        n += 1
        error = np.abs(term / sum_numerator)
        if np.all(error < threshold):
            converged = True

    sum_denominator = np.float64(0)
    n = 1
    converged = False
    while not converged:
        term = (k_radpcm * s_cm * cosd(theta_deg)) ** (2 * n) / math.factorial(n) \
               * np.abs(F_t + 2 ** (n + 2) * R_v0 \
                   / (np.exp((k_radpcm * s_cm * cosd(theta_deg)) ** 2) * cosd(theta_deg))) ** 2 \
               * w_n_isotropic(k_r=2 * k_radpcm * sind(theta_deg),
                               n=n,
                               l_cm=l_cm,
                               alpha=alpha)
        sum_denominator += term.astype(float)
        n += 1
        error = np.abs(term / sum_denominator)
        if np.all(error < threshold):
            converged = True

    S_t = np.abs(F_t) ** 2 * sum_numerator / sum_denominator

    R_ht = R_h + (R_h0 - R_h) * (1 - S_t / S_t0)
    R_vt = R_v + (R_v0 - R_v) * (1 - S_t / S_t0)

    return R_h, R_ht, R_v, R_vt



if validation:
    import matplotlib.pyplot as plt
    plt.style.use('../../custom.mplstyle')

    k_radpcm = 2*np.pi / 3
    s_cm = np.linspace(start=0.001, stop=5/k_radpcm, num=10000)
    epsilon = 50 - 1j * 50
    theta_deg = 80

    R_h = np.zeros_like(s_cm, dtype=complex)
    R_ht = np.zeros_like(s_cm, dtype=complex)
    R_v = np.zeros_like(s_cm, dtype=complex)
    R_vt = np.zeros_like(s_cm, dtype=complex)
    R_h0 = np.zeros_like(s_cm, dtype=complex)
    R_v0 = np.zeros_like(s_cm, dtype=complex)


    for i, s in enumerate(s_cm):
        R_h[i], R_ht[i], R_v[i], R_vt[i] = fresnel_transition(epsilon=epsilon,
                                                              theta_deg=theta_deg,
                                                              k_radpcm=k_radpcm,
                                                              s_cm=s,
                                                              l_cm=4*s,
                                                              alpha=2)
        R_h0[i], R_v0[i] = fresnel(epsilon=epsilon, theta_deg=0)
    fig = plt.figure(num=1, figsize=(10, 5))
    ax = fig.add_subplot(121)
    ax.plot(k_radpcm * s_cm, np.abs(R_h), 'k--', label=r'$R_h(\theta)$')
    ax.plot(k_radpcm * s_cm, np.abs(R_ht), 'b', label=r'$R_{ht}$')
    ax.plot(k_radpcm * s_cm, np.abs(R_h0), 'k-.', label=r'$R_{h0}$')
    ax.set_xlabel(r'$ks$')
    ax.set_ylim([0, 1])
    plt.legend(loc='lower right')

    ax = fig.add_subplot(122)
    ax.plot(k_radpcm * s_cm, np.abs(R_v), 'k--', label=r'$R_v(\theta)$')
    ax.plot(k_radpcm * s_cm, np.abs(R_vt), 'b', label=r'$R_{vt}$')
    ax.plot(k_radpcm * s_cm, np.abs(R_v0), 'k-.', label=r'$R_{v0}$')
    ax.set_xlabel(r'$ks$')
    ax.set_ylim([0, 1])
    plt.legend(loc='lower right')
    plt.suptitle(r'$\varepsilon = $' + f'{epsilon:.0f}' + r', $\theta= $' + f'{theta_deg:.0f}' + f'$^\circ$',
                 fontsize=24)
    plt.tight_layout()
    plt.savefig('../../validation_figures/transition_function.png')
