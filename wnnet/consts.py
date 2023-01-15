import numpy as np
from gslconsts.consts import *

# Conversions for use with the package.

c = GSL_CONST_CGSM_SPEED_OF_LIGHT
k_B = GSL_CONST_CGSM_BOLTZMANN
m_e = GSL_CONST_CGSM_MASS_ELECTRON
u = GSL_CONST_CGSM_UNIFIED_ATOMIC_MASS
N_A = GSL_CONST_NUM_AVOGADRO
hbar = GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR

MeV_to_ergs = 1.0e6 * GSL_CONST_CGSM_ELECTRON_VOLT
ergs_to_MeV = 1.0 / MeV_to_ergs

m_e_in_MeV = m_e * np.power(c, 2) * ergs_to_MeV
m_u_in_MeV = u * np.power(c, 2) * ergs_to_MeV
