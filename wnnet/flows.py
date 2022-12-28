import wnutils.xml as wx
import wnnet.base as fb
import numpy as np


def compute_flows_for_zones(net, zones, nuc_xpath="", reac_xpath=""):
    """A class to compute flows for a set of zones.

    Args:
        ``net``: A wnnet network.

        ``zones`` (:obj:`dict`): A dictionary of `wnutils <https://wnutils.readthedocs.io>`_ *zone data*.

        ``nuc_xpath`` (:obj:`str`, optional): XPath expression
        to select nuclides for flow computations.  Defaults to all
        species.

        ``reac_xpath`` (:obj:`str`, optional): XPath expression
        to select reactions for flow computations.  Defaults to all
        reactions.

    Returns:
        A :obj:`dict` of flows for each zone.  The data for
        each zone are themselves a :obj:`dict` of reactions with each
        item in the dictionary a tuple giving the forward and
        reverse flow.

    """

    nuclides = net.get_nuclides()
    reactions = net.get_reactions()

    dups = net.compute_duplicate_factors(reactions)

    zone_flows = {}

    for zone in zones:
        s_t9 = "t9"
        s_rho = "rho"
        props = zones[zone]["properties"]
        x = zones[zone]["mass fractions"]
        f = {}
        if s_t9 in props and s_rho in props:
            for reaction in net.get_valid_reactions(
                nuc_xpath=nuc_xpath, reac_xpath=reac_xpath
            ):
                my_reaction = reactions[reaction]

                forward, reverse = net.compute_rates_for_reaction(nuclides, my_reaction, float(props[s_t9]), float(props[s_rho]))

                forward *= np.power(
                    float(props[s_rho]), len(my_reaction.nuclide_reactants) - 1
                )
                forward /= dups[reaction][0]
                for sp in my_reaction.nuclide_reactants:
                    tup = net.xml.get_z_a_state_from_nuclide_name(sp)
                    key = (sp, tup[0], tup[1])
                    if key in x:
                        y = x[key] / nuclides[sp]["a"]
                    else:
                        y = 0
                    forward *= y

                if not net.is_weak_reaction(my_reaction):
                    reverse *= np.power(
                        float(props[s_rho]), len(my_reaction.nuclide_products) - 1
                    )
                    reverse /= dups[reaction][1]
                    for sp in my_reaction.nuclide_products:
                        tup = net.xml.get_z_a_state_from_nuclide_name(sp)
                        key = (sp, tup[0], tup[1])
                        if key in x:
                            y = x[key] / nuclides[sp]["a"]
                        else:
                            y = 0
                        reverse *= y
                else:
                    reverse = 0

                f[reaction] = (forward, reverse)

        zone_flows[zone] = f

    return zone_flows
