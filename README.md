# radarscatter
Implements a variety of microwave radar scattering and soil dielectric models for remote sensing applications.

- to install, type ``pip install git+https://github.com/djshiltz/radarscatter``
- to uninstall, type ``pip uninstall radarscatter``

## Dielectric Mixing Models

Computes the complex dielectric constant of a soil-water mixture

### ``'hallikainen_1985'``
Hallikainen et al., 1985, "Microwave Dielectric Behavior of Wet Soil - Part I: Empirical Models and Experimental Observations,"
*IEEE Transactions on Geoscience and Remote Sensing*.

Inputs:
- clay mass fraction [%]
- sand mass fraction [%]
- volumetric moisture [%]
- radar frequency [GHz]

### ``'mironov_2009'``
Mironov et al., 2009, "Physically and Mineralogically Based Spectroscopic Dielectric Model for Moist Soils,"
*IEEE Transactions on Geoscience and Remote Sensing*.

Inputs:
- clay mass fraction [%]
- volumetric moisture [%]
- radar frequency [GHz]

## Radar Scattering Models

Computes the normalized backscattering coefficient of a bare soil surface

### ``'fung_1992'``
Fung et al., 1992, "Backscattering from a Randomly Rough Dielectric Surface," *IEEE Transactions on Geoscience and Remote Sensing*.

The optional transition function for the reflectivity coefficients is from Wu et al., 2001, "A Transition Model for the Reflection Coefficients in Surface Scattering,"
*IEEE Transactions on Geoscience and Remote Sensing*.

The simplified equations for both the IEM model and transition function are listed in Ch. 3 of the textbook: Fung and Chen, 2010, "Microwave
Scattering and Emission Models for Users," *Artech House*.

Inputs:
- radar frequency [GHz]
- radar polarization (hh or vv)
- radar incidence angle [deg]
- complex dielectric constant
- surface RMS height [cm]
- surface correlation length [cm]
- surface autocorrelation function shape parameter (1=exponential, 2=gaussian)
- boolean ``use_transition_function`` flag

### ``'calibrated_iem'``
