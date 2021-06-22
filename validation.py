import numpy as np
import radarscatter
import matplotlib.pyplot as plt
plt.style.use('custom.mplstyle')


## Hallikainen 1985
# Soil properties of Fields 1-5: Table 1, p. 4
sand_pct = [51.51, 41.96, 30.63, 17.16, 5.02]
clay_pct = [13.43, 8.53, 13.48, 19.00, 47.38]
freq_colors = ['red', 'blue', 'green']

# Figure 7, p. 7
fig = plt.figure(figsize=(12, 5))
for field_index, field in enumerate([1, 3, 5]):
    ax = fig.add_subplot(1, 3, field_index + 1)
    for freq_index, freq in enumerate([4.0, 10.0, 18.0]):
        mv_pct = []
        epsilon_real = []
        epsilon_imag = []
        for mv_index, mv in enumerate(np.linspace(0, 50, 1000)):
            ep_r = radarscatter.dielectric(dielectric_model='hallikainen_1985',
                                           freq_GHz=freq,
                                           mv_pct=float(mv),
                                           clay_pct=clay_pct[field - 1],
                                           sand_pct=sand_pct[field - 1])
            epsilon_real.append(np.real(ep_r))
            epsilon_imag.append(-np.imag(ep_r))
            mv_pct.append(mv)
        ax.plot(mv_pct, epsilon_real, color=freq_colors[freq_index], linestyle='solid', label=f'{freq:.0f} GHz')
        ax.plot(mv_pct, epsilon_imag, color=freq_colors[freq_index], linestyle='dashed')
    ax.set_ylim([0, 30])
    ax.set_xticks([0, 10, 20, 30, 40, 50])
    ax.set_xlabel(r'$m_v$ [%]')
    ax.set_ylabel(r'Dielectric Constant $\varepsilon$')
    ax.set_title('Field ' + str(field))
    ax.legend(loc='upper left')
fig.suptitle('Hallikainen et al., 1985 - Fig. 7', fontsize=36)
plt.tight_layout()
plt.savefig('validation_figures/hallikainen_1985.png')


## Mironov 2009 - Reproduce Fig. 7 from Hallikainen et al., 1985 and compare

# Soil properties of Fields 1-5: Table 1, p. 4
sand_pct = [51.51, 41.96, 30.63, 17.16, 5.02]
clay_pct = [13.43, 8.53, 13.48, 19.00, 47.38]
freq_colors = ['red', 'blue', 'green']

# Figure 7, p. 7
fig = plt.figure(figsize=(12, 5))
for field_index, field in enumerate([1, 3, 5]):
    ax = fig.add_subplot(1, 3, field_index + 1)
    for freq_index, freq in enumerate([4.0, 10.0, 18.0]):
        mv_pct = []
        epsilon_real = []
        epsilon_imag = []
        for mv_index, mv in enumerate(np.linspace(0, 50, 1000)):
            ep_r = radarscatter.dielectric(dielectric_model='mironov_2009',
                                           freq_GHz=freq,
                                           mv_pct=float(mv),
                                           clay_pct=clay_pct[field - 1])
            epsilon_real.append(np.real(ep_r))
            epsilon_imag.append(-np.imag(ep_r))
            mv_pct.append(mv)
        ax.plot(mv_pct, epsilon_real, color=freq_colors[freq_index], linestyle='solid', label=f'{freq:.0f} GHz')
        ax.plot(mv_pct, epsilon_imag, color=freq_colors[freq_index], linestyle='dashed')
    ax.set_ylim([0, 30])
    ax.set_xticks([0, 10, 20, 30, 40, 50])
    ax.set_xlabel(r'$m_v$ [%]')
    ax.set_ylabel(r'Dielectric Constant $\varepsilon$')
    ax.set_title('Field ' + str(field))
    ax.legend(loc='upper left')
fig.suptitle('Mironov et al., 2009 - Compare w/ Hallikainen et al., 1985 - Fig. 7', fontsize=36)
plt.tight_layout()
plt.savefig('validation_figures/mironov_2009.png')


## Baghdadi 2016
freq_GHz = 1.0  # arbitrary
lambda_cm = 2.99792458e+8 * 100 / (freq_GHz * 1e+9)
k_radpcm = 2 * np.pi / lambda_cm
ks_array = np.linspace(0.1, 6, 1000)
s_cm_array = ks_array / k_radpcm

fig = plt.figure(figsize=(8, 12), constrained_layout=True)
n_rows = 3
n_cols = 2
for row, pq in enumerate(['hh', 'vv', 'hv']):
    for col, theta_deg in enumerate([25.0, 45.0]):
        index = row * n_cols + col + 1
        ax = fig.add_subplot(n_rows, n_cols, index)
        for mv_pct in [5.0, 15.0, 35.0]:
            sigma_dB_list = []
            for s_cm in s_cm_array:
                sigma_dB_list.append(radarscatter.backscatter(scattering_model='baghdadi_2016',
                                                              freq_GHz=freq_GHz,
                                                              pq=pq,
                                                              theta_deg=theta_deg,
                                                              mv_pct=mv_pct,
                                                              s_cm=float(s_cm)))
            ax.plot(ks_array, sigma_dB_list, label=r'$m_v =$ ' + str(mv_pct) + r' %')
        ax.set_xlabel(r'$ks$')
        ax.set_ylabel(f'$\sigma_{{ {pq} }}$ [dB]')
        ax.set_xlim([0, 6])
        ax.set_ylim([-35, 0])
        ax.set_yticks([-35, -30, -25, -20, -15, -10, -5, 0])
        #ax.legend(loc='lower right')
        ax.set_title(f'$\\theta = {theta_deg}^\circ$')
fig.suptitle('Baghdadi et al., 2016 (Fig. 6)', fontsize=36)
plt.savefig('validation_figures/baghdadi_2016.png')


## Dubois 1995
freq_GHz = 10.0  # arbitrary
lambda_cm = 2.99792458e+8 * 100 / (freq_GHz * 1e+9)
k_radpcm = 2 * np.pi / lambda_cm
ks_array = np.linspace(0.1, 6, 1000)
s_cm_array = ks_array / k_radpcm

fig = plt.figure(figsize=(8, 12), constrained_layout=True)
n_rows = 2
n_cols = 2
for row, pq in enumerate(['hh', 'vv']):
    for col, theta_deg in enumerate([25.0, 45.0]):
        index = row * n_cols + col + 1
        ax = fig.add_subplot(n_rows, n_cols, index)
        for mv_pct in [5.0, 15.0, 35.0]:
            sigma_dB_list = []
            for s_cm in s_cm_array:
                sigma_dB_list.append(radarscatter.backscatter(scattering_model='dubois_1995',
                                                              freq_GHz=freq_GHz,
                                                              pq=pq,
                                                              theta_deg=theta_deg,
                                                              mv_pct=mv_pct,
                                                              s_cm=float(s_cm),
                                                              clay_pct=20.0,
                                                              dielectric_model='mironov_2009'))
            ax.plot(ks_array, sigma_dB_list, label=r'$m_v =$ ' + str(mv_pct) + r' %')
        ax.set_xlabel(r'$ks$')
        ax.set_ylabel(f'$\sigma_{{ {pq} }}$ [dB]')
        ax.set_xlim([0, 6])
        ax.set_ylim([-35, 0])
        ax.set_yticks([-35, -30, -25, -20, -15, -10, -5, 0])
        #ax.legend(loc='lower right')
        ax.set_title(f'$\\theta = {theta_deg}^\circ$')
fig.suptitle('Dubois 1995 (compare w/ Baghdadi 2016)', fontsize=36)
plt.savefig('validation_figures/dubois_1995.png')


## Fung 1992
theta_deg = np.linspace(start=0.01, stop=70, num=100)

# Figure 3.5 (page 58)
fig35 = plt.figure(figsize=(5, 10))
freq_GHz = 5
l_cm = 5
s_cm_list = [0.2, 0.4, 0.8]
linecolors = ['blue', 'black', 'red']
epsilon_real = [6, 36]
epsilon_imag = [0.2, 1.2]
for plot in range(2):
    ax = fig35.add_subplot(211 + plot)
    for j, s_cm in enumerate(s_cm_list):
        sigma_hh = radarscatter.backscatter(scattering_model='fung_1992',
                                            theta_deg=theta_deg,
                                            epsilon_real=epsilon_real[plot],
                                            epsilon_imag=epsilon_imag[plot],
                                            freq_GHz=freq_GHz,
                                            pq='hh',
                                            s_cm=s_cm,
                                            l_cm=l_cm,
                                            alpha=1,
                                            use_transition_function=False)
        sigma_vv = radarscatter.backscatter(scattering_model='fung_1992',
                                            theta_deg=theta_deg,
                                            epsilon_real=epsilon_real[plot],
                                            epsilon_imag=epsilon_imag[plot],
                                            freq_GHz=freq_GHz,
                                            pq='vv',
                                            s_cm=s_cm,
                                            l_cm=l_cm,
                                            alpha=1,
                                            use_transition_function=False)
        ax.plot(theta_deg, sigma_hh, color=linecolors[j], linestyle='solid')
        ax.plot(theta_deg, sigma_vv, color=linecolors[j], linestyle='dashed')
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel(r'$\sigma^0$')
    ax.set_title(r'v & h polarization')
    ax.set_xticks([0, 10, 20, 30, 40, 50, 60, 70])
fig35.suptitle('Fung et al., 1992 (Textbook Fig. 3.5)', fontsize=36)
plt.tight_layout()
plt.savefig('validation_figures/fung_1992_fig3-5.png')

# Fig 3.6 (page 59)
fig36 = plt.figure(figsize=(5, 10))
epsilon_real = 16
epsilon_imag = 1.2
freq_GHz = 5
s_cm = 0.6
l_cm_list = [3, 9, 27, 56]
linecolors = ['red', 'black', 'blue', 'green']
alpha = 1
titles = ['vv Polarization', 'hh Polarization']
pq = ['vv', 'hh']
for plot in range(2):
    ax = fig36.add_subplot(211 + plot)
    for j, l_cm in enumerate(l_cm_list):
        sigma_pp = radarscatter.backscatter(scattering_model='fung_1992',
                                            theta_deg=theta_deg,
                                            epsilon_real=epsilon_real,
                                            epsilon_imag=epsilon_imag,
                                            freq_GHz=freq_GHz,
                                            pq=pq[plot],
                                            s_cm=s_cm,
                                            l_cm=l_cm,
                                            alpha=1,
                                            use_transition_function=False)
        ax.plot(theta_deg, sigma_pp, color=linecolors[j], linestyle='solid')
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel(r'$\sigma^0$')
    ax.set_title(titles[plot])
    ax.set_xticks([0, 10, 20, 30, 40, 50, 60, 70])
fig36.suptitle('Fung et al., 1992 (Textbook Fig. 3.6)', fontsize=36)
plt.tight_layout()
plt.savefig('validation_figures/fung_1992_fig3-6.png')


## Fresnel Equations
from radarscatter.miscellaneous.fresnel_equations import fresnel
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
fig.suptitle('Ulaby et al., 2014 (Fig. 2-19)', fontsize=36)
plt.tight_layout
plt.savefig('validation_figures/fresnel_equations.png')
