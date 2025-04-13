"""This module handles `webnucleo <https://webnucleo.readthedocs.io>`_ nuclear reaction networks."""

import numpy as np
import wnnet.nuc as wn
import wnnet.reac as wr
import wnnet.consts as wc


class Net(wn.Nuc, wr.Reac):
    """A class to store webnucleo networks.

    Args:
        ``file`` (:obj:`str`): A string giving the name of the XML file with\
          the network data.

        ``nuc_xpath`` (:obj:`str`, optional): An XPath expression to select\
         nuclides.  Default is all nuclides.

        ``reac_xpath`` (:obj:`str`, optional):  An XPath expression to select\
        reactions.  Default is all reactions.

    """

    def __init__(self, file, nuc_xpath="", reac_xpath=""):
        wn.Nuc.__init__(self, file, nuc_xpath=nuc_xpath)
        wr.Reac.__init__(self, file, reac_xpath=reac_xpath)
        self.valid_reactions = {}
        self.valid_reactions[("", "")] = self.get_valid_reactions(
            nuc_xpath=nuc_xpath, reac_xpath=reac_xpath
        )

    def compute_q_values(self, nuc_xpath="", reac_xpath=""):
        """A method to compute reaction Q values for valid reactions in the network.

        Args:
            ``nuc_xpath`` (:obj:`str`, optional): An XPath expression to select\
            the nuclides.

            ``reac_xpath`` (:obj:`str`, optional): An XPath expression to\
            select the reactions.

        Returns:
            A :obj:`dict`.  The keys for the dictionary are the reaction\
            strings and the values are the corresponding Q values.

        """

        result = {}

        for reac in self.get_valid_reactions(
            nuc_xpath=nuc_xpath, reac_xpath=reac_xpath
        ):
            tmp = self.compute_reaction_q_value(reac)
            if tmp:
                result[reac] = tmp

        return result

    def get_valid_reactions(self, nuc_xpath="", reac_xpath=""):
        """Method to retrieve the valid reactions in the network.

        Args:
            ``nuc_xpath`` (:obj:`str`, optional):  An XPath expression to\
              select nuclides.  Default is all nuclides.

            ``reac_xpath`` (:obj:`str`, optional):  An XPath expression to\
              select reactions.  Default is all reactions.

        Returns:
            A :obj:`dict` of `wnutils <https://wnutils.readthedocs.io>`_ \
              reactions.

        """

        if (nuc_xpath, reac_xpath) not in self.valid_reactions:
            result = {}
            reactions = self.get_reactions(reac_xpath=reac_xpath)
            for reac in reactions:
                if self.is_valid_reaction(reac, nuc_xpath=nuc_xpath):
                    result[reac] = reactions[reac]
            self.valid_reactions[(nuc_xpath, reac_xpath)] = result

        return self.valid_reactions[(nuc_xpath, reac_xpath)]

    def compute_reaction_q_value(self, name):
        """Method to compute the Q value for a reaction.

        Args:
            ``name`` (:obj:`str`):  A string giving the reaction.

        Returns:
            A :obj:`float` giving the Q value for the reaction or None if the\
             reaction is not valid.

        """

        nuclides = self.get_nuclides()
        reaction = self.get_reactions()[name]
        result = 0
        for spec in reaction.nuclide_reactants:
            if spec not in nuclides:
                return None
            result += nuclides[spec]["mass excess"]
        for spec in reaction.nuclide_products:
            if spec not in nuclides:
                return None
            result -= nuclides[spec]["mass excess"]

        if (
            "positron" in reaction.products
            and "neutrino_e" in reaction.products
        ):
            result -= 2.0 * wc.m_e_in_MeV

        return result

    def _obeys_conservation_laws(self, nuclides, reaction):
        e_p = {"electron": -1, "positron": 1}
        e_lepton = {
            "electron": 1,
            "positron": -1,
            "neutrino_e": 1,
            "anti-neutrino_e": -1,
        }
        mu_lepton = {"neutrino_mu": 1, "anti-neutrino_mu": -1}
        tau_lepton = {"neutrino_tau": 1, "anti-neutrino_tau": -1}

        d_z = 0
        d_a = 0
        d_l_e = 0
        d_l_mu = 0
        d_l_tau = 0

        for spec in reaction.reactants:
            if spec in nuclides:
                d_a += nuclides[spec]["a"]
                d_z += nuclides[spec]["z"]
            if spec in e_p:
                d_z += e_p[spec]
            if spec in e_lepton:
                d_l_e += e_lepton[spec]
            if spec in mu_lepton:
                d_l_mu += mu_lepton[spec]
            if spec in tau_lepton:
                d_l_tau += tau_lepton[spec]
        for spec in reaction.products:
            if spec in nuclides:
                d_a -= nuclides[spec]["a"]
                d_z -= nuclides[spec]["z"]
            if spec in e_p:
                d_z -= e_p[spec]
            if spec in e_lepton:
                d_l_e -= e_lepton[spec]
            if spec in mu_lepton:
                d_l_mu -= mu_lepton[spec]
            if spec in tau_lepton:
                d_l_tau -= tau_lepton[spec]

        return (
            d_a == 0
            and d_z == 0
            and d_l_e == 0
            and d_l_mu == 0
            and d_l_tau == 0
        )

    def is_valid_reaction(self, name, nuc_xpath=""):
        """Method to determine a reaction is valid in the network.

        Args:
            ``name`` (:obj:`str`):  A string giving the reaction.

            ``nuclides`` (:obj:`str`, optional):  An XPath expression\
              to select nuclides.  Default is all nuclides.

        Returns:
            A :obj:`bool` with value True if the reaction is valid and False\
            if not.

        """

        nuclides = self.get_nuclides(nuc_xpath=nuc_xpath)
        reaction = self.get_reactions()[name]

        if not self._obeys_conservation_laws(nuclides, reaction):
            return False

        for spec in reaction.nuclide_reactants + reaction.nuclide_products:
            if spec not in nuclides:
                return False
        return True

    def compute_rates_for_reaction(self, name, t_9, user_funcs=""):
        """Method to compute the forward and reverse rates for a valid reaction.

        Args:
            ``name`` (:obj:`str`):  A string giving the reaction.

            ``t_9`` (:obj:`float`):  The temperature in 10\\ :sup:`9` K at\
             which to compute the rates.

            ``user_funcs`` (:obj:`dict`, optional): A dictionary of\
             user-defined functions associated with a user_rate key.

        Returns:
            A two-element :obj:`tuple` with the first element being the\
            forward rate and the second element being the reverse rate.\
            If the reaction is not valid, returns None.

        """

        if not self.is_valid_reaction(name):
            return None

        reaction = self.get_reactions()[name]
        forward = reaction.compute_rate(t_9, user_funcs=user_funcs)

        if self.is_weak_reaction(name):
            return (forward, 0)

        d_exp = 0

        for spec in reaction.nuclide_reactants:
            d_exp += self.compute_nse_factor(spec, t_9, 1.0)
        for spec in reaction.nuclide_products:
            d_exp -= self.compute_nse_factor(spec, t_9, 1.0)

        if d_exp < -300.0:
            return (forward, 0)

        if d_exp > 300.0:
            return (0, 0)

        tup = self.compute_reaction_duplicate_factors(name)

        return (forward, np.exp(d_exp) * (tup[1] / tup[0]) * forward)

    def compute_rates(self, t_9, nuc_xpath="", reac_xpath="", user_funcs=""):
        """Method to compute the forward and reverse rates for valid\
           reactions in a network.

        Args:
            ``t_9`` (:obj:`float`):  The temperature in 10\\ :sup:`9` K at\
              which to compute the rates.

            ``nuc_xpath`` (:obj:`str`, optional):  An XPath expression\
              to select nuclides.  Default is all nuclides.

            ``reac_xpath`` (:obj:`str`, optional):  An XPath expression\
              to select reactions.  Default is all reactions.

            ``user_funcs`` (:obj:`dict`, optional): A dictionary of\
              user-defined functions associated with a user_rate key.

        Returns:
            A :obj:`dict` containing the rates.  The key is the reaction\
            string while the value is a two-element  :obj:`tuple` with the\
            first element being the forward rate and the second element\
            being the reverse rate.

        """

        v_reactions = self.get_valid_reactions(
            nuc_xpath=nuc_xpath, reac_xpath=reac_xpath
        )

        result = {}

        for reac in v_reactions:
            result[reac] = self.compute_rates_for_reaction(
                reac, t_9, user_funcs=user_funcs
            )

        return result
