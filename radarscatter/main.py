import numpy as np
from radarscatter.scattering_models.baghdadi_2016 import baghdadi_2016
from radarscatter.scattering_models.dubois_1995 import dubois_1995
from radarscatter.scattering_models.fung_1992 import fung_1992
from radarscatter.scattering_models.calibrated_iem import calibrated_iem
from radarscatter.dielectric_models.mironov_2009 import mironov_2009
from radarscatter.dielectric_models.hallikainen_1985 import hallikainen_1985


def dielectric(dielectric_model : str,
               freq_GHz : float,
               mv_pct : float,
               clay_pct : float,
               sand_pct : float = None):
    """
    Computes the complex dielectric constant of a soil-water mixture.

    Required Inputs:
    ----------------

    dielectric_model : float
        Which dielectric mixing model to use, must be in ['hallikainen_1985', 'mironov_2009']

    freq_GHz : float
        Radar frequency [GHz]

    mv_pct : float
        Volumetric soil moisture content (0-100)[%]

    clay_pct : float
        Clay mass fraction (0-100)[%]

    Optional Inputs (depends on model chosen):
    ------------------------------------------

    sand_pct : float
        Sand mass fraction (0-100)[%]

    Output:
    -------

    sigma0_dB
        Normalized radar cross section [dB]
    """

    if dielectric_model not in ['hallikainen_1985', 'mironov_2009']:
        raise ValueError("Model '" + str(dielectric_model) + "' must be in ['hallikainen_1985', 'mironov_2009']")
    if dielectric_model == 'hallikainen_1985':
        if sand_pct is None:
            raise ValueError("Model 'hallikainen_1985' requires 'sand_pct' to be specified")
        return hallikainen_1985(clay_pct=clay_pct,
                                sand_pct=sand_pct,
                                mv_pct=mv_pct,
                                freq_GHz=freq_GHz)
    if dielectric_model == 'mironov_2009':
        return mironov_2009(clay_pct=clay_pct,
                            mv_pct=mv_pct,
                            freq_GHz=freq_GHz)


def backscatter(scattering_model: str,
                pq: str,
                freq_GHz : float,
                theta_deg: float,
                s_cm : float,
                l_cm : float = None,
                alpha : float = None,
                clay_pct: float = None,
                sand_pct: float = None,
                mv_pct: float = None,
                epsilon_real : float = None,
                epsilon_imag : float = None,
                dielectric_model : str = None,
                use_transition_function : bool = None,
                error_threshold : float = 1e-8,
                ):
    """
    Computes normalized radar cross section sigma0_dB using the chosen scattering model

    Required Inputs:
    ----------------

    scattering_model : str
        Which radar scattering model to use, must be in ['dubois_1995', 'baghdadi_2016']

    pq : str
        Radar polarization, must be in ['hh', 'vv', 'hv']

    freq_GHz : float
        Radar frequency [GHz]

    theta_deg : float
        Incidence (sensor zenith) angle [deg]

    s_cm : float
        RMS height of soil surface [cm]

    Optional Inputs (depends on models chosen):
    -------------------------------------------

    l_cm : float
        Correlation length of soil surface [cm]

    alpha : float
        Shape parameter of autocorrelation function [unitless]
            1 - exponential
            2 - gaussian

    clay_pct : float
        Clay mass fraction (0-100)[%]

    sand_pct : float
        Sand mass fraction (0-100)[%]

    mv_pct : float
        Volumetric soil moisture content (0-100)[%]

    epsilon_real : float
        Real part of (relative) dielectric constant [unitless]

    epsilon_imag : float
        Negative imaginary part of (relative) dielectric constant [unitless]
        Always a positive number

    dielectric_model : float
        Which dielectric mixing model to use, must be in ['hallikainen_1985', 'mironov_2009']

    use_transition_function : bool
        Whether to use the transition function outlined in
        Chen et al., 2001 for the IEM

    error_threshold : float
        Error threshold for truncating infinite series in
        the IEM, defaults to 1e-8

    Output:
    -------

    sigma0_dB
        Normalized radar cross section [dB]
    """
    if pq not in ['hh', 'vv', 'hv']:
        raise ValueError("Polarization 'pq' must be in ['hh', 'vv', 'hv']")
    if scattering_model not in ['baghdadi_2016', 'dubois_1995', 'fung_1992', 'calibrated_iem']:
        raise ValueError("Model '" + str(scattering_model) + \
                         "' must be in ['baghdadi_2016', 'dubois_1995', 'fung_1992', 'calibrated_iem']")

    # Common Parameters
    lambda_cm = 2.99792458e+10 / (freq_GHz * 1e+9)
    k_radpcm = 2 * np.pi / lambda_cm

    # Baghdadi 2016
    if scattering_model == 'baghdadi_2016':
        if mv_pct is None:
            raise ValueError("Model 'baghdadi_2016' requires 'mv_pct' to be specified")
        return baghdadi_2016(ks=k_radpcm * s_cm,
                             pq=pq,
                             theta_deg=theta_deg,
                             mv_pct=mv_pct)

    # Dubois 1995
    if scattering_model == 'dubois_1995':
        if epsilon_real is None:
            # user didn't supply dielectric constant, so use dielectric model
            if dielectric_model is None:
                raise ValueError("Please provide a dielectric constant or select a dielectric model")
            epsilon = dielectric(dielectric_model=dielectric_model,
                                 freq_GHz=freq_GHz,
                                 mv_pct=mv_pct,
                                 clay_pct=clay_pct,
                                 sand_pct=sand_pct)
            epsilon_real = float(np.real(epsilon))
        return dubois_1995(ks=k_radpcm * s_cm,
                           lambda_cm=lambda_cm,
                           pq=pq,
                           theta_deg=theta_deg,
                           epsilon_real=epsilon_real)

    # Fung 1992 (original IEM)
    if scattering_model == 'fung_1992':
        if epsilon_real is None or epsilon_imag is None:
            # User did not supply complex dielectric, so use dielectric model
            if dielectric_model is None:
                raise ValueError("Please provide a dielectric constant or select a dielectric model")
            epsilon = dielectric(dielectric_model=dielectric_model,
                                 freq_GHz=freq_GHz,
                                 mv_pct=mv_pct,
                                 clay_pct=clay_pct,
                                 sand_pct=sand_pct)
            epsilon_real = np.real(epsilon)
            epsilon_imag = np.imag(epsilon)
        if l_cm is None:
            raise ValueError("Model 'fung_1992' requires 'l_cm' to be specified")
        if alpha is None:
            raise ValueError("Model 'fung_1992' requires 'alpha' to be specified")
        if use_transition_function is None:
            raise ValueError("Model 'fung_1992' requires 'use_transition_function' to be specified")
        # error threshold is optional, don't check for this
        return fung_1992(theta_deg=theta_deg,
                         epsilon_real=epsilon_real,
                         epsilon_imag=epsilon_imag,
                         freq_GHz=freq_GHz,
                         pq=pq,
                         s_cm=s_cm,
                         l_cm=l_cm,
                         alpha=alpha,
                         use_transition_function=use_transition_function,
                         error_threshold=error_threshold)

    # Calibrated IEM
    if scattering_model == 'calibrated_iem':
        if epsilon_real is None or epsilon_imag is None:
            # User did not supply complex dielectric, so use dielectric model
            epsilon = dielectric(dielectric_model=dielectric_model,
                                 freq_GHz=freq_GHz,
                                 mv_pct=mv_pct,
                                 clay_pct=clay_pct,
                                 sand_pct=sand_pct)
            epsilon_real = np.real(epsilon)
            epsilon_imag = np.imag(epsilon)
        return calibrated_iem(theta_deg=theta_deg,
                              epsilon_real=epsilon_real,
                              epsilon_imag=epsilon_imag,
                              freq_GHz=freq_GHz,
                              pq=pq,
                              s_cm=s_cm,
                              error_threshold=error_threshold)
