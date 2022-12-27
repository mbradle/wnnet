import wnutils.xml as wx
import wnnet.base as fb
import numpy as np


class Flows(fb.Base):
    """A class for computing and graphing flows."""

    def __init__(self, xml):
        self.xml = xml
        self.nuclides = xml.get_nuclide_data()
        self.reactions = xml.get_reaction_data()

    def compute_flows_for_zones(self, zones, induced_nuc_xpath="", induced_reac_xpath=""):
        """A class to compute flows for a set of zones.

        Args:
            ``zones`` (:obj:`dict`): A dictionary of `wnutils <https://wnutils.readthedocs.io>`_ *zone data*.

            ``induced_nuc_xpath`` (:obj:`str`, optional): XPath expression
            to select nuclides for flow computations.  Defaults to all
            species.

            ``induced_reac_xpath`` (:obj:`str`, optional): XPath expression
            to select reactions for flow computations.  Defaults to all
            reactions.

        Returns:
            :obj:`dict`: A dictionary of flows for each zone.  The data for
            each zone are themselvs a :obj:`dict` of reactions with each
            item in the dictionary a tuple giving the forward and
            reverse flow.

        """

        subset_nuclides = {}
        subset_reactions = {}
        if induced_nuc_xpath is None:
            subset_nuclides = self.nuclides
        else:
            subset_nuclides = self.xml.get_nuclide_data(nuc_xpath=induced_nuc_xpath)
        if induced_reac_xpath is None:
            subset_reactions = self.reactions
        else:
            subset_reactions = self.xml.get_reaction_data(reac_xpath=induced_reac_xpath)

        zone_flows = {}

        for zone in zones:
            s_t9 = "t9"
            s_rho = "rho"
            props = zones[zone]["properties"]
            x = zones[zone]["mass fractions"]
            f = {}
            if s_t9 in props and s_rho in props:
                for reaction in self.get_valid_reactions(
                    subset_nuclides, subset_reactions
                ):
                    my_reaction = subset_reactions[reaction]
                    dup_f, dup_r = self.compute_reaction_duplicate_factors(my_reaction)

                    forward = my_reaction.compute_rate(float(props[s_t9]))
                    forward *= np.power(
                        float(props[s_rho]), len(my_reaction.nuclide_reactants) - 1
                    )
                    forward *= dup_f
                    for sp in my_reaction.nuclide_reactants:
                        tup = self.xml.get_z_a_state_from_nuclide_name(sp)
                        key = (sp, tup[0], tup[1])
                        if key in x:
                            y = x[key] / subset_nuclides[sp]["a"]
                        else:
                            y = 0
                        forward *= y

                    if not self.is_weak_reaction(my_reaction):
                        reverse = my_reaction.compute_rate(float(props[s_t9]))
                        reverse *= np.power(
                            float(props[s_rho]), len(my_reaction.nuclide_products) - 1
                        )
                        reverse *= dup_r
                        for sp in my_reaction.nuclide_products:
                            tup = self.xml.get_z_a_state_from_nuclide_name(sp)
                            key = (sp, tup[0], tup[1])
                            if key in x:
                                y = x[key] / subset_nuclides[sp]["a"]
                            else:
                                y = 0
                            reverse *= y
                    else:
                        reverse = 0

                    f[reaction] = (forward, reverse)

            zone_flows[zone] = f

        return zone_flows
