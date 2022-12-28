import wnutils.xml as wx
import numpy as np
from astropy.constants import c, u, N_A, k_B, hbar
import astropy.units as un
from scipy.interpolate import interp1d

class Nuc:
    """A class for handling nuclei and their data."""

    def __init__(self, file):
        self.xml = wx.Xml(file)
        self.nuclides = self.xml.get_nuclide_data()

    def get_nuclides(self, nuc_xpath=""):
        if not nuc_xpath:
            return self.nuclides
        else:
            return self.xml.get_nuclide_data(nuc_xpath=nuc_xpath)

    def compute_nuclear_partition_function(self, nuclide, t9):
        """Method to compute the nuclear partition function for a species.

        Args:
            ``nuclide``: A `wnutils <https://wnutils.readthedocs.io>`_ nuclide.

            ``t9`` (:obj:float): The temperature in 10\ :sup:`9` K at which
            to compute the partition function.

        Returns:
            A :obj:`float` giving the nuclear partition function for the species
            at the input temperature.

        """

        t = nuclide["t9"]
        lg = np.log10(nuclide["partf"])

        if t9 < t[0]:
            return np.power(10.0, lg[0])
        elif t9 > t[len(t) - 1]:
            return np.power(10.0, lg[len(t) - 1])

        if len(t) <= 2:
            f = interp1d(t, lg, kind="linear")
            return np.power(10.0, f(t9))
        else:
            f = interp1d(t, lg, kind="cubic")
            return np.power(10.0, f(t9))

    def compute_quantum_abundance(self, nuclide, t9, rho):
        assert t9 > 0 and rho > 0

        m_u = u * c**2
        m = (m_u.to("MeV").value * nuclide["a"] + nuclide["mass excess"]) * un.MeV

        result = self.compute_nuclear_partition_function(nuclide, t9) / (
            rho * N_A.value
        )

        p1 = m.to("erg").value * k_B.cgs.value * t9 * 1.0e9
        p2 = 2.0 * np.pi * np.power(hbar.cgs.value * c.cgs.value, 2)

        result *= np.power((p1 / p2), 1.5)

        return result

    def compute_binding_energy(self, nuclides, nuclide):
        delta_p = nuclides["h1"]["mass excess"]
        delta_n = nuclides["n"]["mass excess"]

        return (
            nuclide["z"] * delta_p
            + (nuclide["a"] - nuclide["z"]) * delta_n
            - nuclide["mass excess"]
        )

    def compute_NSE_factor(self, nuclides, nuclide, t9, rho):
        B = self.compute_binding_energy(nuclides, nuclide) * un.MeV
        return np.log(self.compute_quantum_abundance(nuclide, t9, rho)) + (
            B.to("erg").value
        ) / ((k_B.cgs * (t9 * 1.0e9 * un.K)).value)

