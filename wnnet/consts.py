"""This module contains the physical constants used in the calculations for
the wnnet package.  webnucleo codes use constants defined in the GNU Scientific
Library (GSL).  For consistency, the wnnet package uses the same constants,
which are defined here.  The constants should be kept up to date with
the definitions in the latest version of GSL."""

import numpy as np

# The GSL cgsm constants are in the file gsl_consts.py.  Those values should be
# kept up to date.

#: Avogadro's number
GSL_CONST_NUM_AVOGADRO = 6.02214199e23

#: Speed of light in vacuum (cm/s)
GSL_CONST_CGS_SPEED_OF_LIGHT = 2.99792458e10

#: Planck's constant divided by 2pi (ergs-s)
GSL_CONST_CGS_PLANCKS_CONSTANT_HBAR = 1.05457162825e-27

#: Electron mass (g)
GSL_CONST_CGS_MASS_ELECTRON = 9.10938188e-28

#: Boltzmann's constant (ergs/K)
GSL_CONST_CGS_BOLTZMANN = 1.3806504e-16

#: Electron volt (ergs)
GSL_CONST_CGS_ELECTRON_VOLT = 1.602176487e-12

#: Amu (g)
GSL_CONST_CGS_UNIFIED_ATOMIC_MASS = 1.660538782e-24

# Conversions for use with the package.

c = GSL_CONST_CGS_SPEED_OF_LIGHT
k_B = GSL_CONST_CGS_BOLTZMANN
m_e = GSL_CONST_CGS_MASS_ELECTRON
u = GSL_CONST_CGS_UNIFIED_ATOMIC_MASS
N_A = GSL_CONST_NUM_AVOGADRO
hbar = GSL_CONST_CGS_PLANCKS_CONSTANT_HBAR

MeV_to_ergs = 1.0e6 * GSL_CONST_CGS_ELECTRON_VOLT
ergs_to_MeV = 1.0 / MeV_to_ergs

m_e_in_MeV = m_e * np.power(c, 2) * ergs_to_MeV
m_u_in_MeV = u * np.power(c, 2) * ergs_to_MeV
