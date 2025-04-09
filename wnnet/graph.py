"""A module to create graphs from network data.

.. _Allowed_keywords:

The routines in this module allow for a number of optional keyword
arguments.  Not all keyword arguments are available for each routine,
however.  The documentation indicates which keyword arguments are
appropriate for a given routine.

The possible keyword arguments are:

.. _induced_nuc_xpath:

    ``induced_nuc_xpath`` (:obj:`str`): An XPath expression to select\
        the subset of nuclides in the graph.  The default is all species\
        in the network.

.. _induced_reac_xpath:

    ``induced_reac_xpath`` (:obj:`str`): An XPath expression to select\
        the subset of reactions in the graph.  The default is all reactions\
        in the network.

.. _flow_type:

    ``flow_type`` (:obj:`str`) A string giving the flow type to be\
        presented.  The possible values are `net`, which shows the forward\
        minus the reverse flow (or the opposite if the reverse flow is larger),\
        and `full`, which shows both the foward and reverse flows.  The\
        default is `net`.

.. _direction:

    ``direction`` (:obj:`str`) A string indicting the reaction directions\
        to show.  Allowed values are `forward`, `reverse`, and `both`\
        (the default, which shows both `forward` and `reverse`).

.. _reaction_color_tuples:

    ``reaction_color_tuples`` (:obj:`list`): A list of :obj:`tuple` objects\
        to select arc colors.  There is a tuple for each reaction type.\
        The first member of the tuple is an XPath expression to select the\
        reaction type while the second member is a string giving the color\
        for that reaction type.  The default is that all arcs are black.

.. _user_funcs:

    ``user_funcs`` (:obj:`dict`): A dictionary of user-defined\
        functions associated with a user_rate key.  The prototype for each\
        user rate function should be (*reaction*, *t_9*), where\
        *t_9* is the temperature in billions of Kelvin and *reaction*\
        is a `wnutils <https://wnutils.readthedocs.io>`_ reaction\
        instance.  Other data can be bound to the function.

.. _zone_user_funcs:

     ``zone_user_funcs`` (:obj:`dict`): A dictionary of user-defined\
        functions associated with a user_rate key.\
        The prototype for each\
        user rate function should be (*reaction*, *t_9*, *zone*), where\
        *t_9* is the temperature in billions of Kelvin and *reaction* and\
        *zone* are `wnutils <https://wnutils.readthedocs.io>`_ reaction and\
        zone instances.  Other data can be bound to the function.

.. _threshold:

    ``threshold`` (:obj:`float`):  The minimum flow (relative to the\
        maximum flow) to be shown on the graph.  The default is 0.01.

.. _scale:
        
    ``scale`` (:obj:`float`):  Scaling factor for the maximum weight arc.\
        The default is 10.

.. _state_scaling:
        
    ``state_scaling`` (:obj:`float`):  Scaling factor for isomeric states.\
        The default is 0.35.

.. _allow_isolated_species:
        
    ``allow_isolated_species`` (:obj:`bool`):  Boolean to choose whether\
        to allow isolated species (ones without incoming or outgoing arcs)\
        in the graph.  The default is *False*.

.. _title_func:

    ``title_func``: A\
       `function <https://docs.python.org/3/library/stdtypes.html#functions>`_\
        that applies the title to the graph.  The function must take\
        one :obj:`float` argument giving the maximum flow.\
        Other data can be bound to the function.  The function must return\
        a :obj:`str` giving the title.  The default is\
        :meth:`wnnet.graph.make_t9_rho_flow_string`.

.. _zone_title_func:
 
    ``zone_title_func``: A\
        `function <https://docs.python.org/3/library/stdtypes.html#functions>`_\
        that applies the title to the graph.  The function must take\
        three arguments.  The first is a :obj:`float` giving the maximum\
        flow.  The second is the zone object corresponding to\
        the graph while the third is the zone label.  Other data can \
        be bound to the function.  The function must return a :obj:`str` \
        giving the title.\
        The default is :meth:`wnnet.graph.make_time_t9_rho_flow_string`\
        for flow graphs and\
        :meth:`wnnet.graph.make_time_t9_rho_current_string` for integrated\
        current graphs.

.. _node_label_func:
        
    ``node_label_func``: A\
      `function  <https://docs.python.org/3/library/stdtypes.html#functions>`_ \
       that applies a label to each node in the graph.  The function\
       must take as argument a species name.  Other data can be bound to\
       the function.  The function must return a :obj:`str`\
       giving the label.\
       The default is :meth:`wnnet.graph.make_node_label`.  The *g_names*\
       argument for the default is bound by the calling routine.

.. _zone_node_label_func:

    ``zone_node_label_func``: A `function \
        <https://docs.python.org/3/library/stdtypes.html#functions>`_ \
        that applies a label to each node in the graph.  The function \
        must take as arguments a species name, a zone, and the zone label. \
        Other data can be bound to \
        the function.  The function must return a :obj:`str` \
        giving the label.  The default is \
        :meth:`wnnet.graph.make_zone_node_label`.

.. _scale_edge_weight_func:
        
    ``scale_edge_weight_func``: A `function \
        <https://docs.python.org/3/library/stdtypes.html#functions>`_ \
        that applies scales each edge weight in the graph.  The function \
        must take four arguments: a dictionary of edge data, the \
        maximum edge weight in the scope of the graph, a scale factor \
        by which to scale the weight (input as `scale`_ to this routine), \
        and a threshold for not including the edge in the graph \
        (input as `threshold`_ to this routine). \
        Other data can be bound to \
        the function.  The function must modify the weight in the \
        edge data and return a :obj:`bool` indicating whether to include \
        the edge in the graph (True) or not (False).\
        The default is :meth:`wnnet.graph.scale_edge_weight`.

.. _graph_attributes:
        
    ``graph_attributes`` (:obj:`dict`):  A dictionary of graphviz attributes\
        for the graph.  The default is *{"outputorder": "edgesfirst"}*.

.. _edge_attributes:
        
    ``edge_attributes`` (:obj:`dict`):  A dictionary of grapvhiz attributes\
        for the edges in the graph.  The default is *{"arrowsize": 0.2}*.

.. _node_attributes:
        
    ``node_attributes`` (:obj:`dict`):  A dictionary of graphviz\
        attributes for the nodes in the graph.  The default is\
        *{ "shape": "box", "fontsize": 16, "style": "filled",\
         "fillcolor": "white" }*.

.. _solar_species:

    ``solar_species`` (:obj:`list`):  A list of species to be\
        considered as the naturally occurring species.  The default is the list\
        returned from :meth:`wnnet.graph.get_solar_species`.

.. _solar_node_attributes:

    ``solar_node_attributes`` (:obj:`dict`):  A dictionary of graphviz\
        attributes to be applied to the solar species in the graph.\
        The default is *{"fillcolor": "yellow", "style": "filled"}*.

.. _special_node_attributes:

    ``special_node_attributes`` (:obj:`dict`):  A dictionary of graphviz\
        attributes to be applied to the special nodes in the graph.\
        The dictionary has as keys the names of the special nodes and as values\
        a dictionary of graphviz properties to be applied to the given special\
        node.

.. _special_edge_attributes:

    ``special_edge_attributes`` (:obj:`dict`):  A dictionary of graphviz\
        attributes to be applied to the special edges in the graph.\
        The dictionary has as keys a :obj:`tuple` *(u, v, reaction)*,\
        where *u* is the source of the directed edge (i.e., arc),\
        *v* is the target,\
        and *reaction* is the reaction from which the arc derives.  The\
        *reaction* tuple element serves as a key in the multi-digraph to\
        distinguish parallel arcs.  The values of the dictionary are each\
        a dictionary of graphviz properties to be applied to the given special\
        edge.

"""

import functools
import networkx as nx
import wnnet.flows as wf
import wnnet.graph_helper as gh


def get_solar_species():
    """A method to return the naturally-occurring solar-system species.

    Returns:
        A :obj:`list`.  A list containing :obj:`str` of the names
        of the naturally-occurring solar-system species.

    """

    return gh.get_solar_species_list()


def make_time_t9_rho_flow_string(f_max, zone, zone_label):
    """The default title function for zone flow graphs.

    Args:
        ``f_max`` (:obj:`float`):  The maximum flow in the scope of the graph.

        ``zone``:  A `wnutils <https://wnutils.readthedocs.io>`_ zone object.

        ``zone_label`` (:obj:`str` or :obj:`tuple`):  The label of the zone.

    Returns:
        A :obj:`str`.  The title string.

    """

    assert zone_label

    props = zone["properties"]
    time = float(props["time"])
    t_9 = float(props["t9"])
    rho = float(props["rho"])

    return f"<time (s) = {gh.fman(time):.2f} x 10<sup>{gh.fexp(time):d}</sup>, \
T<sub>9</sub> = {t_9:.2f}, \
rho (g/cc) = {gh.fman(rho):.2f} x 10<sup>{gh.fexp(rho):d}</sup> (g/cc), \
max. flow = {gh.fman(f_max):.2f} x 10<sup>{gh.fexp(f_max):d}</sup> (s<sup>-1</sup>)>"


def make_time_t9_rho_current_string(f_max, zone, zone_label):
    """The default title function for zone integrated current graphs.

    Args:
        ``f_max`` (:obj:`float`):  The maximum integrated current in the scope\
          of the graph.

        ``zone``:  A `wnutils <https://wnutils.readthedocs.io>`_ zone object.

        ``zone_label`` (:obj:`str` or :obj:`tuple`):  The label of the zone.

    Returns:
        A :obj:`str`.  The title string.

    """

    assert zone_label

    props = zone["properties"]
    time = float(props["time"])
    t_9 = float(props["t9"])
    rho = float(props["rho"])

    return f"<time (s) = {gh.fman(time):.2f} x 10<sup>{gh.fexp(time):d}</sup>, \
T<sub>9</sub> = {t_9:.2f}, \
rho (g/cc) = {gh.fman(rho):.2f} x 10<sup>{gh.fexp(rho):d}</sup> (g/cc), \
max. int. current = {gh.fman(f_max):.2f} x 10<sup>{gh.fexp(f_max):d}</sup>>"


def make_t9_rho_flow_string(f_max, t_9, rho):
    """The default title function for flow graphs.

    Args:
        ``f_max`` (:obj:`float`):  The maximum flow in the scope of the graph.

        ``t_9`` (:obj:`float`):  The temperature in 10\\ :sup:`9` K at which
        the flows are computed.  This value is bound by the calling function.

        ``rho`` (:obj:`float`):  The density in g/cc at which the flows are
        computed.  This value is bound by the calling function.

    Returns:
        A :obj:`str`.  The title string.

    """

    return f"<T<sub>9</sub> = {t_9:.2f}, \
rho (g/cc) = {gh.fman(rho):.2f} x 10<sup>{gh.fexp(rho):d}</sup> (g/cc), \
max. flow = {gh.fman(f_max):.2f} x 10<sup>{gh.fexp(f_max):d}</sup> (s<sup>-1</sup>)>"


def make_node_label(name, g_names):
    """The default node label function.

    Args:
        ``name`` (:obj:`str`):  The species name.

        ``g_names`` (:obj:`dict`):  A dictionary of graphviz names.\
         This dictionary is bound to the function by the calling routine.

    Returns:
        A :obj:`str`.  The node label.

    """
    return g_names[name]


def make_zone_node_label(name, zone, zone_label, g_names):
    """The default zone node label function.

    Args:
        ``name`` (:obj:`str`):  The species name.

        ``zone``:  A `wnutils <https://wnutils.readthedocs.io>`_ zone object.

        ``zone_label`` (:obj:`str` or :obj:`tuple`):  The label of the zone.

        ``g_names`` (:obj:`dict`):  A dictionary of graphviz names.
        This dictionary is bound to the function by the calling routine.

    Returns:
        A :obj:`str`.  The node label.

    """
    assert zone
    assert zone_label

    return make_node_label(name, g_names)


def scale_edge_weight(edge_data, f_max, scale, threshold):
    """The default edge weight scale function.

    Args:
        ``edge_data``:  A dictionary of edge properties.

        ``f_max`` (:obj:`float`):  The maximum edge weight.

        ``scale`` (:obj:`float`):  The factor by which to scale the weight.

        ``threshold`` (:obj:`float`):  The threshold for not including the edge.

    Returns:
        A :obj:`bool` indicating whether to include the edge (True) or not (False).

    """

    if f_max == 0:
        return False

    if "penwidth" in edge_data:
        return True

    keep_edge = True
    _r = edge_data["weight"] / f_max
    if _r >= threshold:
        edge_data["penwidth"] = scale * _r
    else:
        keep_edge = False
    return keep_edge


def create_flow_graph(net, t_9, rho, mass_fractions, **kwargs):
    """A routine to create a flow graph for a given set of mass fractions at
       the input temperature and density.

    Args:
        ``net``: A wnnet network. 

        ``t_9`` (:obj:`float`):  The temperature in 10\\ :sup:`9` K at which to\
         compute the flows. 

        ``rho`` (:obj:`float`):  The density in g/cc at which to compute the\
          flows.
        
        ``mass_fractions`` (:obj:`float`): A\
          `wnutils <https://wnutils.readthedocs.io>`_ dictionary of mass\
          fractions.

        ``**kwargs``: The allowed optional\
          :ref:`keyword <Allowed_keywords>` arguments for this routine are:\
          `induced_nuc_xpath`_, `induced_reac_xpath`_, `flow_type`_,\
          `user_funcs`_, `reaction_color_tuples`_, `threshold`_, `scale`_,\
          `state_scaling`_, `allow_isolated_species`_, `title_func`_,\
          `node_label_func`_, `scale_edge_weight_func`_, `graph_attributes`_,\
          `edge_attributes`_, `node_attributes`_, `solar_species`_,\
          `solar_node_attributes`_, `special_node_attributes`_,\
          `special_edge_attributes`_.


    Returns:
        A `networkx multidigraph \
        <https://networkx.org/documentation/stable/reference/classes/multidigraph.html>`_\
        showing the flows.
     

    """

    my_list = [
        "induced_nuc_xpath",
        "induced_reac_xpath",
        "flow_type",
        "user_funcs",
        "reaction_color_tuples",
        "threshold",
        "scale",
        "state_scaling",
        "allow_isolated_species",
        "title_func",
        "node_label_func",
        "scale_edge_weight_func",
        "graph_attributes",
        "edge_attributes",
        "node_attributes",
        "solar_species",
        "solar_node_attributes",
        "special_node_attributes",
        "special_edge_attributes",
    ]

    my_args = gh.get_keywords(my_list, **kwargs)

    assert my_args["flow_type"] in ("net", "full")

    my_flows = wf.compute_flows(
        net,
        t_9,
        rho,
        mass_fractions,
        reac_xpath=my_args["induced_reac_xpath"],
        user_funcs=my_args["user_funcs"],
    )

    # Get the subset of nuclides to view in the graph.  Get anchors.

    subset_nuclides, anchors = gh.get_subset_and_anchors(
        net.get_nuclides(nuc_xpath=my_args["induced_nuc_xpath"])
    )

    # Title

    if not my_args["title_func"]:
        my_args["title_func"] = functools.partial(
            make_t9_rho_flow_string, t_9=t_9, rho=rho
        )

    # Node label

    if not my_args["node_label_func"]:
        g_names = net.xml.get_graphviz_names(subset_nuclides)
        my_args["node_label_func"] = functools.partial(
            make_node_label, g_names=g_names
        )

    # Scale edge weight

    if not my_args["scale_edge_weight_func"]:
        my_args["scale_edge_weight_func"] = scale_edge_weight

    return gh.create_flow_graph(
        net, my_flows, subset_nuclides, anchors, **my_args
    )


def create_zone_flow_graphs(net, zones, **kwargs):
    """A routine to create flow graphs for a set of zones.

    Args:
        ``net``: A wnnet network. 

        ``zones``: A `wnutils <https://wnutils.readthedocs.io>`_\
         dictionary of zones.

        ``**kwargs``: The allowed optional\
          :ref:`keyword <Allowed_keywords>` arguments for this routine are:\
          `induced_nuc_xpath`_, `induced_reac_xpath`_, `flow_type`_,\
          `zone_user_funcs`_, `reaction_color_tuples`_, `threshold`_,\
          `scale`_, `state_scaling`_, `allow_isolated_species`_,\
          `zone_title_func`_, `zone_node_label_func`_,\
          `scale_edge_weight_func`_,\
          `graph_attributes`_, `edge_attributes`_, `node_attributes`_,\
          `solar_species`_, `solar_node_attributes`_,\
          `special_node_attributes`_, `special_edge_attributes`_.

    Returns:
        A :obj:`dict` of\
        `networkx multidigraphs\
         <https://networkx.org/documentation/stable/reference/classes/multidigraph.html>`_\
         showing the flows.  The keys are the zone labels.
     
    """

    my_list = [
        "induced_nuc_xpath",
        "induced_reac_xpath",
        "flow_type",
        "zone_user_funcs",
        "reaction_color_tuples",
        "threshold",
        "scale",
        "state_scaling",
        "allow_isolated_species",
        "zone_title_func",
        "zone_node_label_func",
        "scale_edge_weight_func",
        "graph_attributes",
        "edge_attributes",
        "node_attributes",
        "solar_species",
        "solar_node_attributes",
        "special_node_attributes",
        "special_edge_attributes",
    ]

    my_args = gh.get_keywords(my_list, **kwargs)

    result = {}

    my_flows = wf.compute_flows_for_zones(
        net,
        zones,
        reac_xpath=my_args["induced_reac_xpath"],
        user_funcs=my_args["zone_user_funcs"],
    )

    subset_nuclides, anchors = gh.get_subset_and_anchors(
        net.get_nuclides(nuc_xpath=my_args["induced_nuc_xpath"])
    )

    g_names = net.xml.get_graphviz_names(subset_nuclides)

    if not my_args["scale_edge_weight_func"]:
        my_args["scale_edge_weight_func"] = scale_edge_weight

    # Loop on zones

    for key, value in zones.items():

        # Title

        if not my_args["zone_title_func"]:
            my_args["title_func"] = functools.partial(
                make_time_t9_rho_flow_string, zone=value, zone_label=key
            )
        else:
            my_args["title_func"] = functools.partial(
                my_args["zone_title_func"], zone=value, zone_label=key
            )

        # Node label

        if not my_args["zone_node_label_func"]:
            my_args["node_label_func"] = functools.partial(
                make_zone_node_label,
                zone=value,
                zone_label=key,
                g_names=g_names,
            )
        else:
            my_args["node_label_func"] = functools.partial(
                my_args["zone_node_label_func"], zone=value, zone_label=key
            )

        # Create graph

        result[key] = gh.create_flow_graph(
            net, my_flows[key], subset_nuclides, anchors, **my_args
        )

    return result


def create_network_graph(net, **kwargs):
    """A routine to create a network graph showing species and reactions\
        among them.

    Args:
        ``net``: A wnnet network. 

        ``**kwargs``: The allowed optional\
          :ref:`keyword <Allowed_keywords>` arguments for this routine are:\
          `induced_nuc_xpath`_, `induced_reac_xpath`_, `direction`_,\
          `reaction_color_tuples`_, `threshold`_, `scale`_, `state_scaling`_,\
          `allow_isolated_species`_, `node_label_func`_, `graph_attributes`_,\
          `edge_attributes`_, `node_attributes`_, `solar_species`_,\
          `solar_node_attributes`_, `special_node_attributes`_,\
          `special_edge_attributes`_.

    Returns:
        A `networkx multidigraph\
        <https://networkx.org/documentation/stable/reference/classes/multidigraph.html>`_\
        showing the network; that is, the nuclides and reactions among them.

    """

    my_list = [
        "induced_nuc_xpath",
        "induced_reac_xpath",
        "direction",
        "reaction_color_tuples",
        "threshold",
        "scale",
        "state_scaling",
        "allow_isolated_species",
        "node_label_func",
        "graph_attributes",
        "node_attributes",
        "edge_attributes",
        "solar_species",
        "solar_node_attributes",
        "special_node_attributes",
        "special_edge_attributes",
    ]

    my_args = gh.get_keywords(my_list, **kwargs)

    assert my_args["direction"] in ("forward", "reverse", "both")

    nuclides = net.get_nuclides()

    # Get the subset of nuclides to view in the graph.  Get anchors.

    val, anchors = gh.get_subset_and_anchors(
        net.get_nuclides(nuc_xpath=my_args["induced_nuc_xpath"])
    )

    d_g = nx.MultiDiGraph()

    # Add nodes.

    for nuc in nuclides:
        d_g.add_node(nuc)

    # Add arcs (reactions).

    gh.add_reactions_to_graph(net, d_g, my_args)

    # Apply attributes

    gh.apply_graph_attributes(d_g, my_args["graph_attributes"])

    gh.apply_node_attributes(d_g, my_args["node_attributes"])

    gh.apply_edge_attributes(d_g, my_args["edge_attributes"])

    gh.apply_solar_node_attributes(
        d_g, my_args["solar_species"], my_args["solar_node_attributes"]
    )

    gh.apply_special_node_attributes(d_g, my_args["special_node_attributes"])

    gh.apply_special_edge_attributes(d_g, my_args["special_edge_attributes"])

    # Remove isolated nodes if desired

    if not my_args["allow_isolated_species"]:
        gh.remove_isolated_nodes(d_g, my_args)

    # Restore anchors

    for anchor in anchors:
        if anchor not in d_g.nodes:
            d_g.add_node(anchor, style="invis")

    # Subgraph

    sub_graph = nx.subgraph(d_g, val)

    # Node label

    if my_args["node_label_func"]:
        _node_label_func = my_args["node_label_func"]
    else:
        g_names = net.xml.get_graphviz_names(list(sub_graph.nodes.keys()))
        _node_label_func = functools.partial(make_node_label, g_names=g_names)

    for node in sub_graph.nodes:
        sub_graph.nodes[node]["pos"] = gh.get_pos(
            net, node, my_args["state_scaling"]
        )
        sub_graph.nodes[node]["label"] = _node_label_func(node)

    gh.color_edges(sub_graph, net, my_args["reaction_color_tuples"])

    return sub_graph


def create_zone_integrated_current_graphs(net, zones, **kwargs):
    """A routine to create integrated currents graphs for a set of zones.

    Args:
        ``net``: A wnnet network. 

        ``zones``: A `wnutils <https://wnutils.readthedocs.io>`_ dictionary\
          of zones.

        ``**kwargs``: The allowed optional\
          :ref:`keyword <Allowed_keywords>` arguments for this routine are:\
          `induced_nuc_xpath`_, `induced_reac_xpath`_,\
          `reaction_color_tuples`_, `threshold`_, `scale`_, `state_scaling`_,\
          `allow_isolated_species`_, `zone_title_func`_,\
          `zone_node_label_func`_,\
          `scale_edge_weight_func`_, `graph_attributes`_, `edge_attributes`_,\
          `node_attributes`_, `solar_species`_, `solar_node_attributes`_,\
          `special_node_attributes`_, `special_edge_attributes`_.

    Returns:
        A :obj:`dict` of `networkx multidigraphs \
        <https://networkx.org/documentation/stable/reference/classes/multidigraph.html>`_ \
        showing the integrated currents.  The keys are the zone labels.

    """

    my_list = [
        "induced_nuc_xpath",
        "induced_reac_xpath",
        "reaction_color_tuples",
        "threshold",
        "scale",
        "state_scaling",
        "allow_isolated_species",
        "zone_title_func",
        "zone_node_label_func",
        "scale_edge_weight_func",
        "graph_attributes",
        "edge_attributes",
        "node_attributes",
        "solar_species",
        "solar_node_attributes",
        "special_node_attributes",
        "special_edge_attributes",
    ]

    my_args = gh.get_keywords(my_list, **kwargs)

    result = {}

    subset_nuclides, anchors = gh.get_subset_and_anchors(
        net.get_nuclides(nuc_xpath=my_args["induced_nuc_xpath"])
    )

    g_names = net.xml.get_graphviz_names(list(net.get_nuclides().keys()))

    if not my_args["scale_edge_weight_func"]:
        my_args["scale_edge_weight_func"] = scale_edge_weight

    for key, value in zones.items():

        # Title

        if not my_args["zone_title_func"]:
            my_args["title_func"] = functools.partial(
                make_time_t9_rho_current_string, zone=value, zone_label=key
            )
        else:
            my_args["title_func"] = functools.partial(
                my_args["zone_title_func"], zone=value, zone_label=key
            )

        # Node label

        if not my_args["zone_node_label_func"]:
            my_args["node_label_func"] = functools.partial(
                make_zone_node_label,
                zone=value,
                zone_label=key,
                g_names=g_names,
            )
        else:
            my_args["node_label_func"] = functools.partial(
                my_args["zone_node_label_func"],
                zone=value,
                zone_label=key,
                g_names=g_names,
            )

        result[key] = gh.create_integrated_current_graph(
            net, value, subset_nuclides, anchors, **my_args
        )

    return result
