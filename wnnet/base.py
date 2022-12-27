import wnutils.xml as wx
import numpy as np
from astropy.constants import c, m_e
from scipy.interpolate import interp1d

class Base:
    """A class for computing and graphing flows."""

    def compute_reaction_Q_value(self, nuclides, reaction):
        result = 0
        for sp in reaction.nuclide_reactants:
            if sp not in nuclides:
                return None
            else:
                result += nuclides[sp]["mass excess"]
        for sp in reaction.nuclide_products:
            if sp not in nuclides:
                return None
            else:
                result -= nuclides[sp]["mass excess"]

        if (
            "positron" in reaction.nuclide_products
            and "neutrino_e" in reaction.nuclide_products
        ):
            E = m_e * c**2
            result -= 2.0 * E.to("MeV").value

        if (
            "positron" in reaction.nuclide_reactants
            and "anti-neutrino_e" in reaction_products
        ):
            E = m_e * c**2
            result += 2.0 * E.to("MeV").value

        return result

    def _compute_duplicate_factor(self, elements):
        dict = {}
        for sp in elements:
            if sp in dict:
                dict[sp] += 1
            else:
                dict[sp] = 1
        result = 1
        for sp in dict:
            result *= np.math.factorial(dict[sp])
        return result

    def compute_reaction_duplicate_factors(self, reaction):
        return (
            self._compute_duplicate_factor(reaction.nuclide_reactants),
            self._compute_duplicate_factor(reaction.nuclide_products),
        )

    def compute_duplicate_factors(self, reactions):
        result = {}
        for r in reactions:
            result[r] = self.compute_reaction_duplicate_factors(reactions[r])
        return result

    def is_weak_reaction(self, reaction):
        result = False
        for sp in reaction.reactants + reaction.products:
            if "electron" in sp or "positron" in sp or "neutrino" in sp:
                result = True
        return result

    def is_valid_reaction(self, nuclides, reaction):
        for sp in reaction.nuclide_reactants + reaction.nuclide_products:
            if sp not in nuclides:
                return False
        return True

    def compute_nuclear_partition_function(self, nuclide, t9):
        "Method to compute the nuclear partition function for a species.

        Args:
            ``nuclide``: A `wnutils <https://wnutils.readthedocs.io>`_ nuclide.

            ``t9`` (:obj:float): The temperature in 10<sup>9</sup> K at which
            to compute the partition function.

        Returns:
            A :obj:float giving the nuclear partition function for the species
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

