import wnutils.xml as wx
import numpy as np
from scipy.interpolate import interp1d


class Reac:
    """A class for handling reactions."""

    def __init__(self, file):
        self.xml = wx.Xml(file)
        self.reactions = self.xml.get_reaction_data()

    def get_reactions(self, reac_xpath=""):
        if not reac_xpath:
            return self.reactions
        else:
            return self.xml.get_reaction_data(reac_xpath=reac_xpath)

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
