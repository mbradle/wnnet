"""This module computes various reaction flows in a network."""

from dataclasses import dataclass
import numpy as np


@dataclass
class _FlowData:
    valid_reactions: dict
    dups: dict
    user_funcs: list
    scale: float = None
    direction: str = None
    order: str = None


def _compute_flows_for_valid_reactions(
    net, t_9, rho, mass_fractions, flow_data: _FlowData
):

    result = {}

    for reaction in flow_data.valid_reactions:
        _reaction = flow_data.valid_reactions[reaction]

        forward, reverse = net.compute_rates_for_reaction(
            reaction, t_9, user_funcs=flow_data.user_funcs
        )

        forward *= np.power(rho, len(_reaction.nuclide_reactants) - 1)
        forward /= flow_data.dups[reaction][0]
        forward *= _compute_abundance_product(
            net, mass_fractions, _reaction.nuclide_reactants
        )

        if not net.is_weak_reaction(reaction):
            reverse *= np.power(rho, len(_reaction.nuclide_products) - 1)
            reverse /= flow_data.dups[reaction][1]
            reverse *= _compute_abundance_product(
                net, mass_fractions, _reaction.nuclide_products
            )
        else:
            reverse = 0

        result[reaction] = (forward, reverse)

    return result


def compute_flows(net, t_9, rho, mass_fractions, **kwargs):
    """A routine to compute flows for a given set of mass fractions at the
       input temperature and density.

    Args:
        ``net``: A wnnet network.

        ``t_9`` (:obj:`float`):  The temperature in 10\\ :sup:`9` K at which
          to compute the flows.

        ``rho`` (:obj:`float`):  The density in g/cc at which to compute the
          flows.

        ``mass_fractions`` (:obj:`float`):
        A `wnutils <https://wnutils.readthedocs.io>`_ dictionary of mass
        fractions.

        ``nuc_xpath`` (:obj:`str`, optional): XPath expression
        to select nuclides for flow computations.  Defaults to all
        species.

        ``reac_xpath`` (:obj:`str`, optional): XPath expression
        to select reactions for flow computations.  Defaults to all
        reactions.

        ``user_funcs`` (:obj:`dict`, optional): A dictionary of user-defined
        functions associated with a user_rate key.
        The prototype for each
        user rate function should be (*reaction*, *t_9*), where
        *t_9* is the temperature in billions of Kelvin and *reaction*
        is a `wnutils <https://wnutils.readthedocs.io>`_ reaction
        instance.  Other data can be bound to the function.

    Returns:
        A :obj:`dict` of reactions with each
        item in the dictionary a tuple giving the forward and
        reverse flow.

    """

    def get_args(**kwargs):
        my_defs = {
            "nuc_xpath": "",
            "reac_xpath": "",
            "user_funcs": "",
        }
        for key in kwargs:
            assert key in my_defs, f"{key} is not an allowed keyword"

        return {**my_defs, **kwargs}

    my_args = get_args(**kwargs)

    valid_reactions = net.get_valid_reactions(
        nuc_xpath=my_args["nuc_xpath"], reac_xpath=my_args["reac_xpath"]
    )

    dups = net.compute_duplicate_factors()

    flow_data = _FlowData(valid_reactions, dups, my_args["user_funcs"])

    return _compute_flows_for_valid_reactions(
        net,
        t_9,
        rho,
        mass_fractions,
        flow_data,
    )


def compute_flows_for_zones(
    net, zones, nuc_xpath="", reac_xpath="", user_funcs=""
):
    """A routine to compute flows for a set of zones.

    Args:
        ``net``: A wnnet network.

        ``zones`` (:obj:`dict`): A dictionary of
         `wnutils <https://wnutils.readthedocs.io>`_ *zone data*.

        ``nuc_xpath`` (:obj:`str`, optional): XPath expression
        to select nuclides for flow computations.  Defaults to all
        species.

        ``reac_xpath`` (:obj:`str`, optional): XPath expression
        to select reactions for flow computations.  Defaults to all
        reactions.

        ``user_funcs`` (:obj:`dict`, optional): A dictionary of user-defined
        functions associated with a user_rate key.
        The prototype for each
        user rate function should be (*reaction*, *t_9*, *zone*), where
        *t_9* is the temperature in billions of Kelvin and *reaction* and
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

    for key, value in zones.items():
        props = value["properties"]
        if "t9" in props and "rho" in props:
            t_9 = float(props["t9"])
            _user_funcs = {}
            if user_funcs:
                for func in user_funcs:
                    _user_funcs[
                        func
                    ] = lambda reaction, t_9, func=func: user_funcs[func](
                        reaction, t_9, value
                    )

            flow_data = _FlowData(valid_reactions, dups, _user_funcs)

            zone_flows[key] = _compute_flows_for_valid_reactions(
                net,
                t_9,
                float(props["rho"]),
                value["mass fractions"],
                flow_data,
            )

    return zone_flows


def _compute_link_flows_for_valid_reactions(
    net,
    t_9,
    rho,
    mass_fractions,
    flow_data: _FlowData,
):
    link_flows = {}

    for reaction in flow_data.valid_reactions:
        tup_array = []
        _reaction = flow_data.valid_reactions[reaction]

        reactants = _reaction.nuclide_reactants
        products = _reaction.nuclide_products

        forward, reverse = net.compute_rates_for_reaction(
            reaction, t_9, user_funcs=flow_data.user_funcs
        )

        if flow_data.direction in ("forward", "both"):
            forward *= np.power(rho, len(_reaction.nuclide_reactants) - 1)
            forward /= flow_data.dups[reaction][0]

            for i, source in enumerate(reactants):
                p_source = _compute_abundance_product(
                    net, mass_fractions, reactants, exclude_index=i
                )
                for target in products:
                    if flow_data.order == "normal":
                        tup = (
                            source,
                            target,
                            forward * p_source * flow_data.scale,
                        )
                    else:
                        tup = (
                            target,
                            source,
                            forward * p_source * flow_data.scale,
                        )
                    tup_array.append(tup)

                if flow_data.direction == "both":
                    for target in reactants:
                        if flow_data.order == "normal":
                            tup = (
                                source,
                                target,
                                -forward * p_source * flow_data.scale,
                            )
                        else:
                            tup = (
                                target,
                                source,
                                -forward * p_source * flow_data.scale,
                            )
                        tup_array.append(tup)

        if not net.is_weak_reaction(reaction):
            if flow_data.direction in ("reverse", "both"):
                reverse *= np.power(rho, len(_reaction.nuclide_products) - 1)
                reverse /= flow_data.dups[reaction][1]

                for i, source in enumerate(products):
                    p_source = _compute_abundance_product(
                        net, mass_fractions, products, exclude_index=i
                    )
                    for target in reactants:
                        if flow_data.order == "normal":
                            tup = (
                                source,
                                target,
                                reverse * p_source * flow_data.scale,
                            )
                        else:
                            tup = (
                                target,
                                source,
                                reverse * p_source * flow_data.scale,
                            )
                        tup_array.append(tup)

                    if flow_data.direction == "both":
                        for target in products:
                            if flow_data.order == "normal":
                                tup = (
                                    source,
                                    target,
                                    -reverse * p_source * flow_data.scale,
                                )
                            else:
                                tup = (
                                    target,
                                    source,
                                    -reverse * p_source * flow_data.scale,
                                )
                            tup_array.append(tup)

        link_flows[reaction] = tup_array

    return link_flows


def compute_link_flows(
    net,
    t_9,
    rho,
    mass_fractions,
    nuc_xpath="",
    reac_xpath="",
    user_funcs="",
    direction="both",
    order="normal",
):
    """A routine to compute link flows for a given set of mass fractions
       at the input temperature and density.

    Args:
        ``net``: A wnnet network.

        ``t_9`` (:obj:`float`):  The temperature in 10\\ :sup:`9` K at which
          to compute the flows.

        ``rho`` (:obj:`float`):  The density in g/cc at which to compute
          the flows.

        ``mass_fractions`` (:obj:`float`):
          A `wnutils <https://wnutils.readthedocs.io>`_ dictionary of mass
          fractions.

        ``nuc_xpath`` (:obj:`str`, optional): XPath expression
        to select nuclides for flow computations.  Defaults to all
        species.

        ``reac_xpath`` (:obj:`str`, optional): XPath expression
        to select reactions for flow computations.  Defaults to all
        reactions.

        ``user_funcs`` (:obj:`dict`, optional): A dictionary of user-defined
        functions associated with a user_rate key.
        The prototype for each
        user rate function should be (*reaction*, *t_9*), where
        *t_9* is the temperature in billions of Kelvin and *reaction*
        is a `wnutils <https://wnutils.readthedocs.io>`_ reaction
        instance.  Other data can be bound to the function.

        ``direction`` (:obj:`str`, optional):  A string indicating the
          direction of the links ("forward", from reactants to products;
          "reverse", from products to reactants; "both", both "forward"
          and "reverse").  Default is "both".

        ``order`` (:obj:`str`, optional):  A string indicating the order
          of the links.  Default is *normal*, in which the *source* and
          *target* of the link are in the time-forward direction of the
          reaction.  For *reversed*, the *source* and *target* are in the
          opposite of the time-forward direction of the reaction such that
          the *target* is the *contribution* to the *source* over some interval.

    Returns:
        A :obj:`dict` of reactions with each
        item in the dictionary an array of three-element :obj:`tuple` objects.
        The tuple elements are the *source*, *target*, and *link flow*.
        The *source* and *target* are determined from the time-forward
        direction of the reaction.

    """

    assert direction in ("forward", "reverse", "both", "all")
    assert order in ("normal", "reversed")

    valid_reactions = net.get_valid_reactions(
        nuc_xpath=nuc_xpath, reac_xpath=reac_xpath
    )

    dups = net.compute_duplicate_factors()

    scale = 1

    flow_data = _FlowData(
        valid_reactions, dups, user_funcs, scale, direction, order
    )

    return _compute_link_flows_for_valid_reactions(
        net,
        t_9,
        rho,
        mass_fractions,
        flow_data,
    )


def compute_link_flows_for_zones(net, zones, **kwargs):
    """A routine to compute link flows for a set of zones.

    Args:
        ``net``: A wnnet network.

        ``zones`` (:obj:`dict`): A dictionary of
          `wnutils <https://wnutils.readthedocs.io>`_ *zone data*.

        ``nuc_xpath`` (:obj:`str`, optional): XPath expression
        to select nuclides for flow computations.  Defaults to all
        species.

        ``reac_xpath`` (:obj:`str`, optional): XPath expression
        to select reactions for flow computations.  Defaults to all
        reactions.

        ``user_funcs`` (:obj:`dict`, optional): A dictionary of user-defined
        functions associated with a user_rate key.
        The prototype for each
        user rate function should be (*reaction*, *t_9*, *zone*), where
        *t_9* is the temperature in billions of Kelvin and *reaction* and
        *zone* are `wnutils <https://wnutils.readthedocs.io>`_ reaction and
        zone instances.  Other data can be bound to the function.

        ``direction`` (:obj:`str`, optional):  A string indicating the
         direction of the links ("forward", from reactants to products;
         "reverse", from products to reactants; "both", both "forward" and
         "reverse").  Default is "both".


        ``include_dt`` (:obj:`bool`, optional):  Boolean determining whether
        to include the *dt* (time interval) in the flow (True) or
        not (False).  Default is False.

        ``order`` (:obj:`str`, optional):  A string indicating the order of
          the links.  Default is *normal*, in which the *source* and *target*
          of the link are in the time-forward direction of the reaction.
          For *reversed*, the *source* and *target* are in the opposite of
          the time-forward direction of the reaction such that the *target*
          is the *contribution* to the *source* over some interval.

    Returns:
        A :obj:`dict` of flow links for each zone.  The data for
        each zone are themselves a :obj:`dict` of reactions with each
        item in the dictionary an array of three-element :obj:`tuple` objects.
        The tuple elements are the *source*, *target*, and *link flow*.
        The *source* and *target* are determined from the time-forward
        direction of the reaction.

    """

    def get_args(**kwargs):
        my_defs = {
            "nuc_xpath": "",
            "reac_xpath": "",
            "user_funcs": "",
            "direction": "both",
            "include_dt": False,
            "order": "normal",
        }
        for key in kwargs:
            assert key in my_defs, f"{key} is not an allowed keyword"

        return {**my_defs, **kwargs}

    my_args = get_args(**kwargs)

    assert my_args["direction"] in ("forward", "reverse", "both")
    assert my_args["order"] in ("normal", "reverse")

    valid_reactions = net.get_valid_reactions(
        nuc_xpath=my_args["nuc_xpath"], reac_xpath=my_args["reac_xpath"]
    )

    dups = net.compute_duplicate_factors()

    zone_link_flows = {}

    for zone in zones:
        props = zones[zone]["properties"]
        if "t9" in props and "rho" in props:
            t_9 = float(props["t9"])
            if my_args["include_dt"]:
                scale = float(props["dt"])
            else:
                scale = 1
            _user_funcs = {}
            if my_args["user_funcs"]:
                for func in my_args["user_funcs"]:
                    _user_funcs[
                        func
                    ] = lambda reaction, t_9, func=func: my_args["user_funcs"][
                        func
                    ](
                        reaction, t_9, zones[zone]
                    )

            flow_data = _FlowData(
                valid_reactions,
                dups,
                _user_funcs,
                scale,
                my_args["direction"],
                my_args["order"],
            )

            zone_link_flows[zone] = _compute_link_flows_for_valid_reactions(
                net,
                t_9,
                float(props["rho"]),
                zones[zone]["mass fractions"],
                flow_data,
            )

    return zone_link_flows


def _compute_abundance_product(net, _x, sp_array, exclude_index=None):
    nuclides = net.get_nuclides()
    result = 1
    for i, _sp in enumerate(sp_array):
        if i is not exclude_index:
            tup = net.xml.get_z_a_state_from_nuclide_name(_sp)
            key = (_sp, tup[0], tup[1])
            if key in _x:
                _y = _x[key] / nuclides[_sp]["a"]
            else:
                _y = 0
            result *= _y
    return result
