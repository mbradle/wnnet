"""This module computes various reaction flows in a network."""

import wnutils.xml as wx
import numpy as np


def _compute_flows_for_valid_reactions(
    net,
    t9,
    rho,
    mass_fractions,
    valid_reactions,
    dups,
    nuc_xpath,
    reac_xpath,
    user_funcs,
):

    result = {}

    for reaction in valid_reactions:
        _reaction = valid_reactions[reaction]

        forward, reverse = net.compute_rates_for_reaction(
            reaction, t9, user_funcs=user_funcs
        )

        forward *= np.power(rho, len(_reaction.nuclide_reactants) - 1)
        forward /= dups[reaction][0]
        forward *= _compute_abundance_product(
            net, mass_fractions, _reaction.nuclide_reactants
        )

        if not net.is_weak_reaction(reaction):
            reverse *= np.power(rho, len(_reaction.nuclide_products) - 1)
            reverse /= dups[reaction][1]
            reverse *= _compute_abundance_product(
                net, mass_fractions, _reaction.nuclide_products
            )
        else:
            reverse = 0

        result[reaction] = (forward, reverse)

    return result


def compute_flows(
    net, t9, rho, mass_fractions, nuc_xpath="", reac_xpath="", user_funcs=""
):
    """A routine to compute flows for a given set of mass fractions at the input temperature and density.

    Args:
        ``net``: A wnnet network.

        ``t9`` (:obj:`float`):  The temperature in 10\ :sup:`9` K at which to compute the flows.

        ``rho`` (:obj:`float`):  The density in g/cc at which to compute the flows.

        ``mass_fractions`` (:obj:`float`): A `wnutils <https://wnutils.readthedocs.io>`_ dictionary of mass fractions.

        ``nuc_xpath`` (:obj:`str`, optional): XPath expression
        to select nuclides for flow computations.  Defaults to all
        species.

        ``reac_xpath`` (:obj:`str`, optional): XPath expression
        to select reactions for flow computations.  Defaults to all
        reactions.

        ``user_funcs`` (:obj:`dict`, optional): A dictionary of user-defined
        functions associated with a user_rate key.
        The prototype for each
        user rate function should be (*reaction*, *t9*), where
        *t9* is the temperature in billions of Kelvin and *reaction*
        is a `wnutils <https://wnutils.readthedocs.io>`_ reaction
        instance.  Other data can be bound to the function.

    Returns:
        A :obj:`dict` of reactions with each
        item in the dictionary a tuple giving the forward and
        reverse flow.

    """

    valid_reactions = net.get_valid_reactions(
        nuc_xpath=nuc_xpath, reac_xpath=reac_xpath
    )

    dups = net.compute_duplicate_factors()

    return _compute_flows_for_valid_reactions(
        net,
        t9,
        rho,
        mass_fractions,
        valid_reactions,
        dups,
        nuc_xpath,
        reac_xpath,
        user_funcs,
    )


def compute_flows_for_zones(
    net, zones, nuc_xpath="", reac_xpath="", user_funcs=""
):
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

        ``user_funcs`` (:obj:`dict`, optional): A dictionary of user-defined
        functions associated with a user_rate key.
        The prototype for each
        user rate function should be (*reaction*, *t9*, *zone*), where
        *t9* is the temperature in billions of Kelvin and *reaction* and
        *zone* are `wnutils <https://wnutils.readthedocs.io>`_ reaction and
        zone instances.  Other data can be bound to the function.


    Returns:
        A :obj:`dict` of flows for each zone.  The data for
        each zone are themselves a :obj:`dict` of reactions with each
        item in the dictionary a tuple giving the forward and
        reverse flow.

    """

    zone_flows = {}

    valid_reactions = net.get_valid_reactions(
        nuc_xpath=nuc_xpath, reac_xpath=reac_xpath
    )

    dups = net.compute_duplicate_factors()

    for zone in zones:
        s_t9 = "t9"
        s_rho = "rho"
        _zone = zones[zone]
        props = _zone["properties"]
        if s_t9 in props and s_rho in props:
            t9 = float(props[s_t9])
            rho = float(props[s_rho])
            _user_funcs = {}
            if user_funcs:
                for func in user_funcs:
                    _user_funcs[
                        func
                    ] = lambda reaction, t9, func=func: user_funcs[func](
                        reaction, t9, _zone
                    )
            zone_flows[zone] = _compute_flows_for_valid_reactions(
                net,
                t9,
                rho,
                _zone["mass fractions"],
                valid_reactions,
                dups,
                nuc_xpath,
                reac_xpath,
                _user_funcs,
            )

    return zone_flows


def _compute_link_flows_for_valid_reactions(
    net,
    t9,
    rho,
    mass_fractions,
    valid_reactions,
    dups,
    scale,
    nuc_xpath,
    reac_xpath,
    user_funcs,
    direction,
    order,
):
    nuclides = net.get_nuclides()

    link_flows = {}

    for reaction in valid_reactions:
        tup_array = []
        _reaction = valid_reactions[reaction]

        reactants = _reaction.nuclide_reactants
        products = _reaction.nuclide_products

        forward, reverse = net.compute_rates_for_reaction(
            reaction, t9, user_funcs=user_funcs
        )

        if direction == "forward" or direction == "both":
            forward *= np.power(rho, len(_reaction.nuclide_reactants) - 1)
            forward /= dups[reaction][0]

            for i in range(len(reactants)):
                source = reactants[i]
                p_source = _compute_abundance_product(
                    net, mass_fractions, reactants, exclude_index=i
                )
                for target in products:
                    if order == "normal":
                        tup = (source, target, forward * p_source * scale)
                    else:
                        tup = (target, source, forward * p_source * scale)
                    tup_array.append(tup)

                if direction == "both":
                    for target in reactants:
                        if order == "normal":
                            tup = (source, target, -forward * p_source * scale)
                        else:
                            tup = (target, source, -forward * p_source * scale)
                        tup_array.append(tup)

        if not net.is_weak_reaction(reaction):
            if direction == "reverse" or direction == "both":
                reverse *= np.power(rho, len(_reaction.nuclide_products) - 1)
                reverse /= dups[reaction][1]

                for i in range(len(products)):
                    source = products[i]
                    p_source = _compute_abundance_product(
                        net, mass_fractions, products, exclude_index=i
                    )
                    for target in reactants:
                        if order == "normal":
                            tup = (source, target, reverse * p_source * scale)
                        else:
                            tup = (target, source, reverse * p_source * scale)
                        tup_array.append(tup)

                    if direction == "both":
                        for target in products:
                            if order == "normal":
                                tup = (
                                    source,
                                    target,
                                    -reverse * p_source * scale,
                                )
                            else:
                                tup = (
                                    target,
                                    source,
                                    -reverse * p_source * scale,
                                )
                            tup_array.append(tup)

        link_flows[reaction] = tup_array

    return link_flows


def compute_link_flows(
    net,
    t9,
    rho,
    mass_fractions,
    nuc_xpath="",
    reac_xpath="",
    user_funcs="",
    direction="both",
    order="normal",
):
    """A routine to compute link flows for a given set of mass fractions at the input temperature and density.

    Args:
        ``net``: A wnnet network.

        ``t9`` (:obj:`float`):  The temperature in 10\ :sup:`9` K at which to compute the flows.

        ``rho`` (:obj:`float`):  The density in g/cc at which to compute the flows.

        ``mass_fractions`` (:obj:`float`): A `wnutils <https://wnutils.readthedocs.io>`_ dictionary of mass fractions.

        ``nuc_xpath`` (:obj:`str`, optional): XPath expression
        to select nuclides for flow computations.  Defaults to all
        species.

        ``reac_xpath`` (:obj:`str`, optional): XPath expression
        to select reactions for flow computations.  Defaults to all
        reactions.

        ``user_funcs`` (:obj:`dict`, optional): A dictionary of user-defined
        functions associated with a user_rate key.
        The prototype for each
        user rate function should be (*reaction*, *t9*), where
        *t9* is the temperature in billions of Kelvin and *reaction*
        is a `wnutils <https://wnutils.readthedocs.io>`_ reaction
        instance.  Other data can be bound to the function.

        ``direction`` (:obj:`str`, optional):  A string indicating the direction of the links ("forward", from reactants to products; "reverse", from products to reactants; "both", both "forward" and "reverse").  Default is "both".

        ``order`` (:obj:`str`, optional):  A string indicating the order of the links.  Default is *normal*, in which the *source* and *target* of the link are in the time-forward direction of the reaction.  For *reversed*, the *source* and *target* are in the opposite of the time-forward direction of the reaction such that the *target* is the *contribution* to the *source* over some interval.

    Returns:
        A :obj:`dict` of reactions with each
        item in the dictionary an array of three-element :obj:`tuple` objects.
        The tuple elements are the *source*, *target*, and *link flow*.
        The *source* and *target* are determined from the time-forward
        direction of the reaction.

    """

    assert (
        direction == "forward"
        or direction == "reverse"
        or direction == "both"
        or direction == "all"
    )
    assert order == "normal" or order == "reversed"

    valid_reactions = net.get_valid_reactions(
        nuc_xpath=nuc_xpath, reac_xpath=reac_xpath
    )

    dups = net.compute_duplicate_factors()

    scale = 1

    return _compute_link_flows_for_valid_reactions(
        net,
        t9,
        rho,
        mass_fractions,
        valid_reactions,
        dups,
        scale,
        nuc_xpath,
        reac_xpath,
        user_funcs,
        direction,
        order,
    )


def compute_link_flows_for_zones(
    net,
    zones,
    nuc_xpath="",
    reac_xpath="",
    user_funcs="",
    direction="both",
    include_dt=False,
    order="normal",
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

        ``user_funcs`` (:obj:`dict`, optional): A dictionary of user-defined
        functions associated with a user_rate key.
        The prototype for each
        user rate function should be (*reaction*, *t9*, *zone*), where
        *t9* is the temperature in billions of Kelvin and *reaction* and
        *zone* are `wnutils <https://wnutils.readthedocs.io>`_ reaction and
        zone instances.  Other data can be bound to the function.

        ``direction`` (:obj:`str`, optional):  A string indicating the direction of the links ("forward", from reactants to products; "reverse", from products to reactants; "both", both "forward" and "reverse").  Default is "both".

        ``include_dt`` (:obj:`bool`, optional):  Boolean determining whether to include the *dt* (time interval) in the flow (True) or not (False).  Default is False.

        ``order`` (:obj:`str`, optional):  A string indicating the order of the links.  Default is *normal*, in which the *source* and *target* of the link are in the time-forward direction of the reaction.  For *reversed*, the *source* and *target* are in the opposite of the time-forward direction of the reaction such that the *target* is the *contribution* to the *source* over some interval.

    Returns:
        A :obj:`dict` of flow links for each zone.  The data for
        each zone are themselves a :obj:`dict` of reactions with each
        item in the dictionary an array of three-element :obj:`tuple` objects.
        The tuple elements are the *source*, *target*, and *link flow*.
        The *source* and *target* are determined from the time-forward
        direction of the reaction.

    """

    assert (
        direction == "forward" or direction == "reverse" or direction == "both"
    )
    assert order == "normal" or order == "reversed"

    nuclides = net.get_nuclides()
    reactions = net.get_reactions()

    s_t9 = "t9"
    s_rho = "rho"
    s_dt = "dt"

    valid_reactions = net.get_valid_reactions(
        nuc_xpath=nuc_xpath, reac_xpath=reac_xpath
    )

    dups = net.compute_duplicate_factors()

    zone_link_flows = {}

    for zone in zones:
        props = zones[zone]["properties"]
        f = {}
        if s_t9 in props and s_rho in props:
            if include_dt:
                scale = float(props[s_dt])
            else:
                scale = 1
            _user_funcs = {}
            if user_funcs:
                for func in user_funcs:
                    _user_funcs[
                        func
                    ] = lambda reaction, t9, func=func: user_funcs[func](
                        reaction, t9, zones[zone]
                    )
            f = _compute_link_flows_for_valid_reactions(
                net,
                float(props[s_t9]),
                float(props[s_rho]),
                zones[zone]["mass fractions"],
                valid_reactions,
                dups,
                scale,
                nuc_xpath,
                reac_xpath,
                _user_funcs,
                direction,
                order,
            )

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
