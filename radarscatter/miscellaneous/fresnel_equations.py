import numpy as np
import radarscatter
from radarscatter.miscellaneous.trig import sind, cosd


validation = False

def fresnel(epsilon, theta_deg):
    """
    Computes complex reflection coefficients R_h and R_v for an interface using
    the equations from the textbook "Microwave Radar and Radiometric Remote Sensing"
    by Ulaby and Moore, page 62.

    Note that in most contexts, R_h and R_v refer to the squared magnitude of these
    quantities.  For the IEM, we want the complex valued version, which is normally
    indicated by a lower case 'r'.

    In the text, these quantities are labeled rho_h and rho_v.  I simplified the
    equations by plugging the expression for cosine of theta_2 into the equations
    for rho, and assuming that medium 1 is air (e1 = 1).

    The sign of R_v is flipped because when Fung does his simplified forms, he assumes
    the Fresnel equations take that form.  If I don't include this, the simplified expressions
    used in his textbook don't give the correct result.

    Inputs:
    -------
    epsilon_real, epsilon_imag : float
        Electrical permittivity ratio (soil to air)

    theta_deg : float
        Incidence (sensor zenith) angle [deg]

    Outputs:
    --------
    R_h, R_v
        Reflection coefficients in h and v polarization
    """

    R_h =  (cosd(theta_deg) - np.sqrt(epsilon - sind(theta_deg)**2)) \
         / (cosd(theta_deg) + np.sqrt(epsilon - sind(theta_deg)**2))

    R_v = -(np.sqrt(epsilon - sind(theta_deg)**2) - epsilon * cosd(theta_deg)) \
        /  (np.sqrt(epsilon - sind(theta_deg)**2) + epsilon * cosd(theta_deg))

    return R_h, R_v


if validation:
    import matplotlib.pyplot as plt
    plt.style.use('../../custom.mplstyle')

    epsilon_real = 50
    epsilon_imag = [0, 25, 50]
    colors = ['green', 'red', 'cyan']

    fig = plt.figure()
    ax = fig.add_subplot()
    for i, epsilon_imag in enumerate(epsilon_imag):
        theta_deg = np.linspace(start=0,
                                stop=90,
                                num=1000,
                                endpoint=True)
        R_h, R_v = fresnel(epsilon=epsilon_real - 1j * epsilon_imag,
                           theta_deg=theta_deg)
        ax.plot(theta_deg, np.abs(R_h), color=colors[i], linestyle='dashed',
                label=r'$|R_h|$ ($\varepsilon = $' + f'{epsilon_real:.0f} - j{epsilon_imag:.0f})')
        ax.plot(theta_deg, np.abs(R_v), color=colors[i], linestyle='solid',
                label=r'$|R_v|$ ($\varepsilon = $' + f'{epsilon_real:.0f} - j{epsilon_imag:.0f})')
    ax.set_xlabel(r'Incidence Angle $\theta$ [deg]')
    plt.legend()
    plt.savefig('../../validation_figures/fresnel_equations.png')
