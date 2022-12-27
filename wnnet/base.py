import wnutils.xml as wx
import numpy as np
from astropy.constants import c, m_e


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

    def compute_duplicate_factor(self, elements):
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
            self.compute_duplicate_factor(reaction.nuclide_reactants),
            self.compute_duplicate_factor(reaction.nuclide_products),
        )

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
