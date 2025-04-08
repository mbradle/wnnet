"""This module handles `webnucleo <https://webnucleo.readthedocs.io>`_ collections of reactions."""

import math
import wnutils.xml as wx


class Reac:
    """A class for handling reactions and their data.

    Args:
        ``file`` (:obj:`str`): A string giving the XML file name with the\
           reaction data.

        ``reac_xpath`` (:obj:`str`, optional):  An XPath expression to\
         select reactions.  Default is all reactions.

    """

    def __init__(self, file, reac_xpath=""):
        self.xml = wx.Xml(file)
        self.reactions = {}
        self.reactions[""] = self.xml.get_reaction_data(reac_xpath=reac_xpath)

    def get_reactions(self, reac_xpath=""):
        """Method to return a collection of reactions.

        Args:
            ``reac_xpath`` (:obj:`str`, optional): An XPath expression\
             to select the reactions.  Default is all reactions in the\
             underlying collection of reactions.

        Returns:
            A :obj:`dict` containing `wnutils <https://wnutils.readthedocs.io>`_ reactions.

        """

        if reac_xpath not in self.reactions:
            self.reactions[reac_xpath] = {}

            new_reactions = self.xml.get_reaction_data(reac_xpath=reac_xpath)

            for key, value in new_reactions.items():
                if key in self.reactions[""]:
                    self.reactions[reac_xpath][key] = value

        return self.reactions[reac_xpath]

    def _compute_duplicate_factor(self, elements):
        my_dict = {}
        for species in elements:
            if species in my_dict:
                my_dict[species] += 1
            else:
                my_dict[species] = 1
        result = 1
        for val in my_dict.values():
            result *= math.factorial(val)
        return result

    def compute_reaction_duplicate_factors(self, name):
        """Method to compute the duplicate reaction element factors for a reaction.

        Args:
            ``name`` (:obj:`str`): A string giving the reaction.

        Returns:
            A two-element :obj:`tuple`.  The first element is the duplicate\
            factor for the forward reaction.  The second element is the\
            duplicate factor for the reverse reaction.

        """

        reaction = self.get_reactions()[name]
        return (
            self._compute_duplicate_factor(reaction.nuclide_reactants),
            self._compute_duplicate_factor(reaction.nuclide_products),
        )

    def compute_duplicate_factors(self, reac_xpath=""):
        """Method to compute the duplicate factors for reactions in the reaction collection.

        Args:
            ``reac_xpath`` (:obj:`str`, optional): An XPath expression\
            to select reactions.  Default is all reactions.

        Returns:
            A :obj:`dict` containing the factors.  The keys are reaction\
            strings.  The values are two-element :obj:`tuple` objects\
            with the first element the duplicate factor for the\
            forward reaction and the second element the duplicate factor\
            for the reverse reaction.

        """

        result = {}
        for reac in self.get_reactions(reac_xpath=reac_xpath):
            result[reac] = self.compute_reaction_duplicate_factors(reac)
        return result

    def is_weak_reaction(self, name):
        """Method to determine if a reaction is a weak reaction or not.

        Args:
            ``name`` (:obj:`str`): A string giving the reaction.

        Returns:
            A :obj:`bool` with value True if the reaction is a weak reaction\
            and False if not.

        """

        reaction = self.get_reactions()[name]
        result = False
        for _sp in reaction.reactants + reaction.products:
            if "electron" in _sp or "positron" in _sp or "neutrino" in _sp:
                result = True
        return result
