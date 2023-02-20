"""This module constants for use in the wnnet package.  The module
also includes all constants from
`gslconsts <https://gslconsts.readthedocs.io/en/latest/gslconsts.html#module-gslconsts.consts>`_.
"""

import numpy as np
from gslconsts.consts import *

# Conversions for use with the package.

#: Speed of light in vacuum (cm/s)
c = GSL_CONST_CGSM_SPEED_OF_LIGHT

#: Boltzmann's constant (ergs/K)
k_B = GSL_CONST_CGSM_BOLTZMANN

#: Electron mass (g)
m_e = GSL_CONST_CGSM_MASS_ELECTRON

#: Amu mass (g)
u = GSL_CONST_CGSM_UNIFIED_ATOMIC_MASS

#: Avogadro's number
N_A = GSL_CONST_NUM_AVOGADRO

#: Planck's constant divided by 2pi (ergs/s)
hbar = GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR

#: Conversion factor from MeVs to ergs
MeV_to_ergs = 1.0e6 * GSL_CONST_CGSM_ELECTRON_VOLT

#: Conversion factor from ergs to MeVs
ergs_to_MeV = 1.0 / MeV_to_ergs

#: Electron rest mass (MeV)
m_e_in_MeV = m_e * np.power(c, 2) * ergs_to_MeV

#: Amu rest mass (MeV)
m_u_in_MeV = u * np.power(c, 2) * ergs_to_MeV
