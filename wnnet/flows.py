"""This module computes various reaction flows in a network."""

import wnutils.xml as wx
import numpy as np


def compute_flows_for_zones(net, zones, nuc_xpath="", reac_xpath=""):
    """A routine to compute flows for a set of zones.

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

    valid_reactions = net.get_valid_reactions(
        nuc_xpath=nuc_xpath, reac_xpath=reac_xpath
    )

    dups = net.compute_duplicate_factors()

    zone_flows = {}

    for zone in zones:
        s_t9 = "t9"
        s_rho = "rho"
        props = zones[zone]["properties"]
        x = zones[zone]["mass fractions"]
        f = {}
        if s_t9 in props and s_rho in props:
            for reaction in valid_reactions:
                my_reaction = reactions[reaction]

                forward, reverse = net.compute_rates_for_reaction(
                    reaction, float(props[s_t9])
                )

                forward *= np.power(
                    float(props[s_rho]), len(my_reaction.nuclide_reactants) - 1
                )
                forward /= dups[reaction][0]
                forward *= _compute_abundance_product(
                    net, x, my_reaction.nuclide_reactants
                )

                if not net.is_weak_reaction(reaction):
                    reverse *= np.power(
                        float(props[s_rho]), len(my_reaction.nuclide_products) - 1
                    )
                    reverse /= dups[reaction][1]
                    reverse *= _compute_abundance_product(
                        net, x, my_reaction.nuclide_products
                    )
                else:
                    reverse = 0

                f[reaction] = (forward, reverse)

        zone_flows[zone] = f

    return zone_flows


def compute_link_flows_for_zones(
    net, zones, nuc_xpath="", reac_xpath="", direction="both", include_dt=False
):
    """A routine to compute link flows for a set of zones.

    Args:
        ``net``: A wnnet network.

        ``zones`` (:obj:`dict`): A dictionary of `wnutils <https://wnutils.readthedocs.io>`_ *zone data*.

        ``nuc_xpath`` (:obj:`str`, optional): XPath expression
        to select nuclides for flow computations.  Defaults to all
        species.

        ``reac_xpath`` (:obj:`str`, optional): XPath expression
        to select reactions for flow computations.  Defaults to all
        reactions.

        ``direction`` (:obj:`str`, optional):  A string indicating the direction of the links ("forward", from reactants to products; "reverse", from products to reactants; "both", both "forward" and "reverse").  Default is "both".

        ``include_dt`` (:obj:`bool`, optional):  Boolean determining whether to include the *dt* (time interval) in the flow (True) or not (False).  Default is False.

    Returns:
        A :obj:`dict` of flow links for each zone.  The data for
        each zone are themselves a :obj:`dict` of reactions with each
        item in the dictionary an array of three-element :obj:`tuple` objects.
        The tuple elements are the *source*, *target*, and *link flow*.
        The *source* and *target* are determined from the time-forward
        direction of the reaction.

    """

    nuclides = net.get_nuclides()
    reactions = net.get_reactions()

    valid_reactions = net.get_valid_reactions(
        nuc_xpath=nuc_xpath, reac_xpath=reac_xpath
    )

    dups = net.compute_duplicate_factors()

    zone_link_flows = {}

    for zone in zones:
        s_t9 = "t9"
        s_rho = "rho"
        s_dt = "dt"
        props = zones[zone]["properties"]
        x = zones[zone]["mass fractions"]
        f = {}
        if s_t9 in props and s_rho in props:
            for reaction in valid_reactions:
                tup_array = []
                my_reaction = reactions[reaction]

                reactants = my_reaction.nuclide_reactants
                products = my_reaction.nuclide_products

                forward, reverse = net.compute_rates_for_reaction(
                    reaction, float(props[s_t9])
                )

                if direction == "forward" or direction == "both":
                    forward *= np.power(
                        float(props[s_rho]), len(my_reaction.nuclide_reactants) - 1
                    )
                    forward /= dups[reaction][0]
                    if include_dt:
                        forward *= float(props[s_dt])

                    for i in range(len(reactants)):
                        source = reactants[i]
                        p_source = _compute_abundance_product(
                            net, x, reactants, exclude_index=i
                        )
                        for target in products:
                            tup = (source, target, forward * p_source)
                            tup_array.append(tup)

                if not net.is_weak_reaction(reaction):
                    if direction == "reverse" or direction == "both":
                        reverse *= np.power(
                            float(props[s_rho]), len(my_reaction.nuclide_products) - 1
                        )
                        reverse /= dups[reaction][1]
                        if include_dt:
                            reverse *= dt

                        for i in range(len(products)):
                            source = products[i]
                            p_source = _compute_abundance_product(
                                net, x, products, exclude_index=i
                            )
                            for target in reactants:
                                tup = (source, target, reverse * p_source)
                                tup_array.append(tup)

                f[reaction] = tup_array

        zone_link_flows[zone] = f

    return zone_link_flows


def _compute_abundance_product(net, x, sp_array, exclude_index=None):
    nuclides = net.get_nuclides()
    result = 1
    for i in range(len(sp_array)):
        if i is not exclude_index:
            tup = net.xml.get_z_a_state_from_nuclide_name(sp_array[i])
            key = (sp_array[i], tup[0], tup[1])
            if key in x:
                y = x[key] / nuclides[sp_array[i]]["a"]
            else:
                y = 0
            result *= y
    return result
