"""A module to create graphs from network data.

.. _Allowed_keywords:

The routines in this module allow for a number of optional keyword
arguments.  Not all keyword arguments are available for each routine,
however.  The doucmentation indicates which keyword arguments are
appropriate for a given routine.

The possible keyword arguments are:

.. _induced_nuc_xpath:

    * ``induced_nuc_xpath`` (:obj:`str`): An XPath expression to select\
        the subset of nuclides in the graph.  The default is all species\
        in the network.

.. _induced_reac_xpath:

    * ``induced_reac_xpath`` (:obj:`str`): An XPath expression to select\
        the subset of reactions in the graph.  The default is all reactions\
        in the network.

.. _flow_type:

    * ``flow_type`` (:obj:`str`) A string giving the flow type to be\
        presented.  The possible values are `net`, which shows the forward\
        minus the reverse flow (or the opposite if the reverse flow is larger),\
        and `full`, which shows both the foward and reverse flows.  The\
        default is `net`.

.. _direction:

    * ``direction`` (:obj:`str`) A string indicting the reaction directions\
        to show.  Allowed values are `forward`, `reverse`, and `both`\
        (the default, which shows both `forward` and `reverse`).

.. _reaction_color_tuples:

    * ``reaction_color_tuples`` (:obj:`list`): A list of :obj:`tuple` objects\
        to select arc colors.  There is a tuple for each reaction type.\
        The first member of the tuple is an XPath\
        expression to select the reaction type while the second member is a\
        string giving the color for that reaction type.  The default is that\
        all arcs are black.

.. _user_funcs:

    * ``user_funcs`` (:obj:`dict`): A dictionary of user-defined\
        functions associated with a user_rate key.  The prototype for each\
        user rate function should be (*reaction*, *t_9*), where\
        *t_9* is the temperature in billions of Kelvin and *reaction*\
        is a `wnutils <https://wnutils.readthedocs.io>`_ reaction\
        instance.  Other data can be bound to the function.

.. _zone_user_funcs:

    *  ``zone_user_funcs`` (:obj:`dict`): A dictionary of user-defined\
        functions associated with a user_rate key.\
        The prototype for each\
        user rate function should be (*reaction*, *t_9*, *zone*), where\
        *t_9* is the temperature in billions of Kelvin and *reaction* and\
        *zone* are `wnutils <https://wnutils.readthedocs.io>`_ reaction and\
        zone instances.  Other data can be bound to the function.\

.. _threshold:

    * ``threshold`` (:obj:`float`):  The minimum flow (relative to the\
        maximum flow) to be shown on the graph.  The default is 0.01.

.. _scale:
        
    * ``scale`` (:obj:`float`):  Scaling factor for the maximum weight arc.\
        The default is 10.

.. _state_scaling:
        
    * ``state_scaling`` (:obj:`float`):  Scaling factor for isomeric states.\
        The default is 0.35.

.. _allow_isolated_species:
        
    * ``allow_isolated_species`` (:obj:`bool`):  Boolean to choose whether\
        to allow isolated species (ones without incoming or outgoing arcs)\
        in the graph.  The default is *False*.

.. _title_func:

    * ``title_func``: A\
       `function <https://docs.python.org/3/library/stdtypes.html#functions>`_\
        that applies the title to the graph.  The function must take\
        one :obj:`float` argument giving the maximum flow.\
        Other data can be bound to the function.  The function must return\
        a :obj:`str` giving the title.  The default is\
        :meth:`wnnet.graph.make_t9_rho_flow_string`.

.. _zone_title_func:
 
    *   ``zone_title_func``: A `function \
        <https://docs.python.org/3/library/stdtypes.html#functions>`_ \
        that applies the title to the graph.  The function must take \
        three arguments.  The first is a :obj:`float` giving the maximum\
        flow.  The second is the zone object corresponding to\
        the graph while the third is the zone label.  Other data can \
        be bound to the function.  The function must return a :obj:`str` \
        giving the title.  \
        The default is :meth:`wnnet.graph.make_time_t9_rho_flow_string`
        for flow graphs and.
        :meth:`wnnet.graph.make_time_t9_rho_current_string` for integrated
        current graphs.

.. _node_label_func:
        
    * ``node_label_func``: A `function \
        <https://docs.python.org/3/library/stdtypes.html#functions>`_ \
        that applies a label to each node in the graph.  The function \
        must take as argument a species name.  Other data can be bound to \
        the function.  The function must return a :obj:`str` \
        giving the label.  \
        The default is :meth:`wnnet.graph.make_node_label`.

.. _zone_node_label_func:

    *  ``zone_node_label_func``: A `function \
        <https://docs.python.org/3/library/stdtypes.html#functions>`_ \
        that applies a label to each node in the graph.  The function \
        must take as arguments a species name, a zone, and the zone label. \
        Other data can be bound to \
        the function.  The function must return a :obj:`str` \
        giving the label.  The default is \
        :meth:`wnnet.graph.make_zone_node_label`.

.. _scale_edge_weight_func:
        
    *  ``scale_edge_weight_func``: A `function \
        <https://docs.python.org/3/library/stdtypes.html#functions>`_ \
        that applies scales each edge weight in the graph.  The function \
        must take as four arguments: a dictionary of edge data, the \
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
        
    * ``graph_attributes`` (:obj:`dict`):  A dictionary of graphviz attributes\
        for the graph.  The default is *{"outputorder": "edgesfirst"}*.

.. _edge_attributes:
        
    * ``edge_attributes`` (:obj:`dict`):  A dictionary of grapvhiz attributes\
        for the edges in the graph.  The default is *{"arrowsize": 0.2}*.

.. _node_attributes:
        
    * ``node_attributes`` (:obj:`dict`):  A dictionary of graphviz\
        attributes for the nodes in the graph.  The default is\
        *{ "shape": "box", "fontsize": 16, "style": "filled",\
         "fillcolor": "white" }*.

.. _solar_species:

    * ``solar_species`` (:obj:`list`):  A list of species to be\
        considered as the naturally occurring species.  The default is the list\
        returned from :meth:`wnnet.graph.get_solar_species`.

.. _solar_node_attributes:

    * ``solar_node_attributes`` (:obj:`dict`):  A dictionary of graphviz\
        attributes to be applied to the solar species in the graph.\
        The default is *{"fillcolor": "yellow", "style": "filled"}*.

.. _special_node_attributes:

    * ``special_node_attributes`` (:obj:`dict`):  A dictionary of graphviz\
        attributes to be applied to the special nodes in the graph.\
        The dictionary has as keys the names of the special nodes and as values\
        a dictionary of graphviz properties to be applied to the given special\
        node.

"""

import operator
from math import floor, log10
import functools
import networkx as nx
import wnnet.flows as wf


def _get_keywords(kw_list, **kwargs):

    for key in kwargs:
        assert key in kw_list

    def_kwargs = {}

    def_kwargs["induced_nuc_xpath"] = ""
    def_kwargs["induced_reac_xpath"] = ""
    def_kwargs["flow_type"] = "net"
    def_kwargs["direction"] = "both"
    def_kwargs["reaction_color_tuples"] = None
    def_kwargs["user_funcs"] = None
    def_kwargs["zone_user_funcs"] = None
    def_kwargs["threshold"] = 0.01
    def_kwargs["scale"] = 10
    def_kwargs["state_scaling"] = 0.35
    def_kwargs["allow_isolated_species"] = False
    def_kwargs["title_func"] = None
    def_kwargs["zone_title_func"] = None
    def_kwargs["node_label_func"] = None
    def_kwargs["zone_node_label_func"] = None
    def_kwargs["scale_edge_weight_func"] = None
    def_kwargs["graph_attributes"] = {"outputorder": "edgesfirst"}
    def_kwargs["node_attributes"] = {
        "shape": "box",
        "fontsize": 16,
        "style": "filled",
        "fillcolor": "white",
    }
    def_kwargs["edge_attributes"] = {"arrowsize": 0.2}
    def_kwargs["solar_species"] = get_solar_species()
    def_kwargs["solar_node_attributes"] = {
        "fillcolor": "yellow",
        "style": "filled",
    }
    def_kwargs["special_node_attributes"] = None

    result = {}

    for k_w in kw_list:
        if k_w in kwargs:
            if isinstance(def_kwargs[k_w], dict):
                result[k_w] = {**def_kwargs[k_w], **kwargs[k_w]}
            else:
                result[k_w] = kwargs[k_w]
        else:
            result[k_w] = def_kwargs[k_w]

    return result


def get_solar_species():
    """A method to return the naturally-occurring solar-system species.

    Returns:
        A :obj:`list`.  A list containing :obj:`str` of the names
        of the naturally-occurring solar-system species.

    """

    return [
        "h1",
        "h2",
        "he3",
        "he4",
        "li6",
        "li7",
        "be9",
        "b10",
        "b11",
        "c12",
        "c13",
        "n14",
        "n15",
        "o16",
        "o17",
        "o18",
        "f19",
        "ne20",
        "ne21",
        "ne22",
        "na23",
        "mg24",
        "mg25",
        "mg26",
        "al27",
        "si28",
        "si29",
        "si30",
        "p31",
        "s32",
        "s33",
        "s34",
        "s36",
        "cl35",
        "cl37",
        "ar36",
        "ar38",
        "ar40",
        "k39",
        "k40",
        "k41",
        "ca40",
        "ca42",
        "ca43",
        "ca44",
        "ca46",
        "ca48",
        "sc45",
        "ti46",
        "ti47",
        "ti48",
        "ti49",
        "ti50",
        "v50",
        "v51",
        "cr50",
        "cr52",
        "cr53",
        "cr54",
        "mn55",
        "fe54",
        "fe56",
        "fe57",
        "fe58",
        "co59",
        "ni58",
        "ni60",
        "ni61",
        "ni62",
        "ni64",
        "cu63",
        "cu65",
        "zn64",
        "zn66",
        "zn67",
        "zn68",
        "zn70",
        "ga69",
        "ga71",
        "ge70",
        "ge72",
        "ge73",
        "ge74",
        "ge76",
        "as75",
        "se74",
        "se76",
        "se77",
        "se78",
        "se80",
        "se82",
        "br79",
        "br81",
        "kr78",
        "kr80",
        "kr82",
        "kr83",
        "kr84",
        "kr86",
        "rb85",
        "rb87",
        "sr84",
        "sr86",
        "sr87",
        "sr88",
        "y89",
        "zr90",
        "zr91",
        "zr92",
        "zr94",
        "zr96",
        "nb93",
        "mo92",
        "mo94",
        "mo95",
        "mo96",
        "mo97",
        "mo98",
        "mo100",
        "ru96",
        "ru98",
        "ru99",
        "ru100",
        "ru101",
        "ru102",
        "ru104",
        "rh103",
        "pd102",
        "pd104",
        "pd105",
        "pd106",
        "pd108",
        "pd110",
        "ag107",
        "ag109",
        "cd106",
        "cd108",
        "cd110",
        "cd111",
        "cd112",
        "cd113",
        "cd114",
        "cd116",
        "in113",
        "in115",
        "sn112",
        "sn114",
        "sn115",
        "sn116",
        "sn117",
        "sn118",
        "sn119",
        "sn120",
        "sn122",
        "sn124",
        "sb121",
        "sb123",
        "te120",
        "te122",
        "te123",
        "te124",
        "te125",
        "te126",
        "te128",
        "te130",
        "i127",
        "xe124",
        "xe126",
        "xe128",
        "xe129",
        "xe130",
        "xe131",
        "xe132",
        "xe134",
        "xe136",
        "cs133",
        "ba130",
        "ba132",
        "ba134",
        "ba135",
        "ba136",
        "ba137",
        "ba138",
        "la138",
        "la139",
        "ce136",
        "ce138",
        "ce140",
        "ce142",
        "pr141",
        "nd142",
        "nd143",
        "nd144",
        "nd145",
        "nd146",
        "nd148",
        "nd150",
        "sm144",
        "sm147",
        "sm148",
        "sm149",
        "sm150",
        "sm152",
        "sm154",
        "eu151",
        "eu153",
        "gd152",
        "gd154",
        "gd155",
        "gd156",
        "gd157",
        "gd158",
        "gd160",
        "tb159",
        "dy156",
        "dy158",
        "dy160",
        "dy161",
        "dy162",
        "dy163",
        "dy164",
        "ho165",
        "er162",
        "er164",
        "er166",
        "er167",
        "er168",
        "er170",
        "tm169",
        "yb168",
        "yb170",
        "yb171",
        "yb172",
        "yb173",
        "yb174",
        "yb176",
        "lu175",
        "lu176",
        "hf174",
        "hf176",
        "hf177",
        "hf178",
        "hf179",
        "hf180",
        "ta180",
        "ta181",
        "w180",
        "w182",
        "w183",
        "w184",
        "w186",
        "re185",
        "re187",
        "os184",
        "os186",
        "os187",
        "os188",
        "os189",
        "os190",
        "os192",
        "ir191",
        "ir193",
        "pt190",
        "pt192",
        "pt194",
        "pt195",
        "pt196",
        "pt198",
        "au197",
        "hg196",
        "hg198",
        "hg199",
        "hg200",
        "hg201",
        "hg202",
        "hg204",
        "tl203",
        "tl205",
        "pb204",
        "pb206",
        "pb207",
        "pb208",
        "bi209",
        "th232",
        "u235",
        "u238",
    ]


def _fexp(_x):
    # Add 0.01 for rounding for :.2f mantissa formating
    return int(floor(log10(abs(_x)) + 0.01)) if _x != 0.0 else 0


def _fman(_x):
    return _x / 10.0 ** _fexp(_x)


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

    return f"<time (s) = {_fman(time):.2f} x 10<sup>{_fexp(time):d}</sup>, \
T<sub>9</sub> = {t_9:.2f}, \
rho (g/cc) = {_fman(rho):.2f} x 10<sup>{_fexp(rho):d}</sup> (g/cc), \
max. flow = {_fman(f_max):.2f} x 10<sup>{_fexp(f_max):d}</sup> (s<sup>-1</sup>)>"


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

    return f"<time (s) = {_fman(time):.2f} x 10<sup>{_fexp(time):d}</sup>, \
T<sub>9</sub> = {t_9:.2f}, \
rho (g/cc) = {_fman(rho):.2f} x 10<sup>{_fexp(rho):d}</sup> (g/cc), \
max. int. current = {_fman(f_max):.2f} x 10<sup>{_fexp(f_max):d}</sup>>"


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
rho (g/cc) = {_fman(rho):.2f} x 10<sup>{_fexp(rho):d}</sup> (g/cc), \
max. flow = {_fman(f_max):.2f} x 10<sup>{_fexp(f_max):d}</sup> (s<sup>-1</sup>)>"


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

    keep_edge = True
    _r = edge_data["weight"] / f_max
    if _r >= threshold:
        edge_data["penwidth"] = scale * _r
    else:
        keep_edge = False
    return keep_edge


def _color_edges(my_graph, net, color_tuples):
    color = {}
    for reaction in net.get_reactions():
        color[reaction] = "black"

    if color_tuples:
        for color_tup in color_tuples:
            for reaction in net.get_reactions(reac_xpath=color_tup[0]):
                color[reaction] = color_tup[1]

        for edge in my_graph.edges:
            my_graph.edges[edge]["color"] = color[
                my_graph.edges[edge]["reaction"]
            ]


def _get_pos(coll, name, state_scaling):
    _z, _a, state = coll.xml.get_z_a_state_from_nuclide_name(name)
    _n = _a - _z
    if state == "g":
        _z -= state_scaling
    elif state == "m":
        _z += state_scaling
    return (float(_n), float(_z))


def _get_subset_and_anchors(nuclides):
    species_subset = []
    my_dict = {}
    z_dict = {}
    for key, value in nuclides.items():
        species_subset.append(key)
        my_dict[(value["z"], value["a"])] = key
        if value["z"] not in z_dict:
            z_dict[value["z"]] = [value["a"]]
        else:
            z_dict[value["z"]].append(value["a"])

    z_array = []

    for key, value in z_dict.items():
        z_array.append(key)
        value.sort()

    z_array.sort()

    anchors = []
    anchors.append(my_dict[(z_array[0], z_dict[z_array[0]][0])])
    anchors.append(my_dict[(z_array[0], z_dict[z_array[0]][-1])])
    anchors.append(my_dict[(z_array[-1], z_dict[z_array[-1]][0])])
    anchors.append(my_dict[(z_array[-1], z_dict[z_array[-1]][-1])])

    return (species_subset, anchors)


def _apply_graph_attributes(my_graph, graph_attributes):
    if graph_attributes:
        for key, value in graph_attributes.items():
            my_graph.graph[key] = value


def _apply_node_attributes(my_graph, node_attributes):
    if node_attributes:
        for key, value in node_attributes.items():
            for node in my_graph.nodes:
                my_graph.nodes[node][key] = value


def _apply_edge_attributes(my_graph, edge_attributes):
    if edge_attributes:
        for key, value in edge_attributes.items():
            for edge in my_graph.edges:
                my_graph.edges[edge][key] = value


def _apply_solar_node_attributes(
    my_graph, solar_species, solar_node_attributes
):

    for node in solar_species:
        if node in my_graph.nodes:
            for key, value in solar_node_attributes.items():
                my_graph.nodes[node][key] = value


def _apply_special_node_attributes(my_graph, special_node_attributes):
    if special_node_attributes:
        for node in special_node_attributes:
            for key, value in special_node_attributes[node].items():
                my_graph.nodes[node][key] = value


def _create_flow_graph(net, my_flows, subset_nuclides, anchors, **my_args):

    reactions = net.get_reactions()

    d_g = nx.MultiDiGraph()

    for nuc in net.get_nuclides():
        d_g.add_node(nuc)

    for key, value in my_flows.items():
        tup = value

        if my_args["flow_type"] == "full":
            if tup[0] > 0:
                for reactant in reactions[key].nuclide_reactants:
                    for product in reactions[key].nuclide_products:
                        d_g.add_edge(
                            reactant, product, weight=tup[0], reaction=key
                        )

            if tup[1] > 0:
                for product in reactions[key].nuclide_products:
                    for reactant in reactions[key].nuclide_reactants:
                        d_g.add_edge(
                            product, reactant, weight=tup[1], reaction=key
                        )

        elif my_args["flow_type"] == "net":
            net_flow = tup[0] - tup[1]

            if net_flow > 0:
                for reactant in reactions[key].nuclide_reactants:
                    for product in reactions[key].nuclide_products:
                        d_g.add_edge(
                            reactant, product, weight=net_flow, reaction=key
                        )

            if net_flow < 0:
                for product in reactions[key].nuclide_products:
                    for reactant in reactions[key].nuclide_reactants:
                        d_g.add_edge(
                            product, reactant, weight=-net_flow, reaction=key
                        )

    # Apply attributes

    _apply_graph_attributes(d_g, my_args["graph_attributes"])

    _apply_node_attributes(d_g, my_args["node_attributes"])

    _apply_edge_attributes(d_g, my_args["edge_attributes"])

    _apply_solar_node_attributes(
        d_g, my_args["solar_species"], my_args["solar_node_attributes"]
    )

    _apply_special_node_attributes(d_g, my_args["special_node_attributes"])

    # Subgraph and maximum flow within subgraph

    sub_graph = nx.subgraph(d_g, subset_nuclides)

    my_weights = nx.get_edge_attributes(sub_graph, "weight")

    # Set penwidth.  Remove edges that are below threshold

    f_max = 0

    if len(my_weights) > 0:
        f_max = max(my_weights.items(), key=operator.itemgetter(1))[1]

        if not my_args["scale_edge_weight_func"]:
            _scale_edge_weight_func = functools.partial(
                scale_edge_weight,
                f_max=f_max,
                scale=my_args["scale"],
                threshold=my_args["threshold"],
            )
        else:
            _scale_edge_weight_func = functools.partial(
                my_args["scale_edge_weight_func"],
                f_max=f_max,
                scale=my_args["scale"],
                threshold=my_args["threshold"],
            )

        remove_edges = []
        for edge in d_g.edges:
            if not _scale_edge_weight_func(d_g.get_edge_data(*edge)):
                remove_edges.append(edge)

        d_g.remove_edges_from(remove_edges)

    # Remove isolated nodes if desired

    if not my_args["allow_isolated_species"]:
        isolated_nodes = list(nx.isolates(d_g))
        for node in isolated_nodes:
            if (
                node not in my_args["solar_species"]
                and node not in my_args["special_node_attributes"]
            ):
                d_g.remove_node(node)

    # Restore anchors

    for anchor in anchors:
        if anchor not in d_g.nodes:
            d_g.add_node(anchor, style="invis")

    # Get new subset

    sub_graph_2 = nx.subgraph(d_g, subset_nuclides)

    for node in sub_graph_2.nodes:
        sub_graph_2.nodes[node]["pos"] = _get_pos(
            net, node, my_args["state_scaling"]
        )
        sub_graph_2.nodes[node]["label"] = my_args["node_label_func"](node)

    # Title

    d_g.graph["label"] = my_args["title_func"](f_max)

    _color_edges(sub_graph_2, net, my_args["reaction_color_tuples"])

    return sub_graph_2


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
          `solar_node_attributes`_, `special_node_attributes`_.


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
    ]

    my_args = _get_keywords(my_list, **kwargs)

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

    subset_nuclides, anchors = _get_subset_and_anchors(
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

    return _create_flow_graph(
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
          `zone_title_func`_, `zone_node_label_func`_, `scale_edge_weight_func`_,\
          `graph_attributes`_, `edge_attributes`_, `node_attributes`_,\
          `solar_species`_, `solar_node_attributes`_, `special_node_attributes`_.

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
    ]

    my_args = _get_keywords(my_list, **kwargs)

    result = {}

    my_flows = wf.compute_flows_for_zones(
        net,
        zones,
        reac_xpath=my_args["induced_reac_xpath"],
        user_funcs=my_args["zone_user_funcs"],
    )

    subset_nuclides, anchors = _get_subset_and_anchors(
        net.get_nuclides(nuc_xpath=my_args["induced_nuc_xpath"])
    )

    g_names = net.xml.get_graphviz_names(subset_nuclides)

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

        result[key] = _create_flow_graph(
            net, my_flows[key], subset_nuclides, anchors, **my_args
        )

    return result


def create_nuclides_graph(nuc, **kwargs):
    """A routine to create a graph showing species.

    Args:
        ``nuc``: A wnnet nuclide object.

        ``**kwargs``: The allowed optional\
          :ref:`keyword <Allowed_keywords>` arguments for this routine are:\
          `induced_nuc_xpath`_, `state_scaling`_, `node_label_func`_,\
          `graph_attributes`_, `node_attributes`_, `solar_species`_,\
          `solar_node_attributes`_, and `special_node_attributes`_.

    Returns:
        A \
       `networkx multidigraph \
        <https://networkx.org/documentation/stable/reference/classes/multidigraph.html>`_\
        showing the nuclides.
     
    """

    my_list = [
        "induced_nuc_xpath",
        "state_scaling",
        "node_label_func",
        "graph_attributes",
        "node_attributes",
        "solar_species",
        "solar_node_attributes",
        "special_node_attributes",
    ]

    my_args = _get_keywords(my_list, **kwargs)

    d_g = nx.MultiDiGraph()

    for species in nuc.get_nuclides(nuc_xpath=my_args["induced_nuc_xpath"]):
        d_g.add_node(species)

    # Apply attributes

    _apply_graph_attributes(d_g, my_args["graph_attributes"])

    _apply_node_attributes(d_g, my_args["node_attributes"])

    _apply_solar_node_attributes(
        d_g, my_args["solar_species"], my_args["solar_node_attributes"]
    )

    _apply_special_node_attributes(d_g, my_args["special_node_attributes"])

    # Node label

    if my_args["node_label_func"]:
        _node_label_func = my_args["node_label_func"]
    else:
        g_names = nuc.xml.get_graphviz_names(list(d_g.nodes.keys()))
        _node_label_func = functools.partial(make_node_label, g_names=g_names)

    for node in d_g.nodes:
        d_g.nodes[node]["pos"] = _get_pos(nuc, node, my_args["state_scaling"])
        d_g.nodes[node]["label"] = _node_label_func(node)

    return d_g


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
          `solar_node_attributes`_, `special_node_attributes`_.

    Returns:
        A `networkx multidigraph \
        <https://networkx.org/documentation/stable/reference/classes/multidigraph.html>`_ \
        `showing the network; that is, the nuclides and reactions among them.
     
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
    ]

    my_args = _get_keywords(my_list, **kwargs)

    assert my_args["direction"] in ("forward", "reverse", "both")

    nuclides = net.get_nuclides()
    reactions = net.get_reactions(reac_xpath=my_args["induced_reac_xpath"])

    # Get the subset of nuclides to view in the graph.  Get anchors.

    val, anchors = _get_subset_and_anchors(
        net.get_nuclides(nuc_xpath=my_args["induced_nuc_xpath"])
    )

    d_g = nx.MultiDiGraph()

    for nuc in nuclides:
        d_g.add_node(nuc)

    for key, value in reactions.items():
        if my_args["direction"] in ("forward", "both"):
            for reactant in value.nuclide_reactants:
                for product in value.nuclide_products:
                    d_g.add_edge(reactant, product, reaction=key)

        if not net.is_weak_reaction(key) and (
            my_args["direction"] in ("reverse", "both")
        ):
            for product in value.nuclide_products:
                for reactant in value.nuclide_reactants:
                    d_g.add_edge(product, reactant, reaction=key)

    # Apply attributes

    _apply_graph_attributes(d_g, my_args["graph_attributes"])

    _apply_node_attributes(d_g, my_args["node_attributes"])

    _apply_edge_attributes(d_g, my_args["edge_attributes"])

    _apply_solar_node_attributes(
        d_g, my_args["solar_species"], my_args["solar_node_attributes"]
    )

    _apply_special_node_attributes(d_g, my_args["special_node_attributes"])

    # Remove isolated nodes if desired

    if not my_args["allow_isolated_species"]:
        d_g.remove_nodes_from(list(nx.isolates(d_g)))

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
        sub_graph.nodes[node]["pos"] = _get_pos(
            net, node, my_args["state_scaling"]
        )
        sub_graph.nodes[node]["label"] = _node_label_func(node)

    _color_edges(sub_graph, net, my_args["reaction_color_tuples"])

    return sub_graph


def _create_integrated_current_graph(
    net, zone, subset_nuclides, anchors, **my_args
):
    reactions = net.get_reactions()

    props = zone["properties"]

    d_g = nx.MultiDiGraph()

    for nuc in net.get_nuclides():
        d_g.add_node(nuc)

    f_currents = {}

    for prop in props:
        if isinstance(prop, tuple):
            if prop[0] == "flow current":
                f_currents[prop[1]] = float(props[prop])

    for key, value in f_currents.items():

        if value > 0:
            for reactant in reactions[key].nuclide_reactants:
                for product in reactions[key].nuclide_products:
                    d_g.add_edge(reactant, product, weight=value, reaction=key)

        if value < 0:
            for product in reactions[key].nuclide_products:
                for reactant in reactions[key].nuclide_reactants:
                    d_g.add_edge(
                        product, reactant, weight=-value, reaction=key
                    )

    # Apply attributes

    _apply_graph_attributes(d_g, my_args["graph_attributes"])

    _apply_node_attributes(d_g, my_args["node_attributes"])

    _apply_edge_attributes(d_g, my_args["edge_attributes"])

    _apply_solar_node_attributes(
        d_g, my_args["solar_species"], my_args["solar_node_attributes"]
    )

    _apply_special_node_attributes(d_g, my_args["special_node_attributes"])

    # Subgraph and maximum flow within subgraph

    sub_graph = nx.subgraph(d_g, subset_nuclides)

    my_weights = nx.get_edge_attributes(sub_graph, "weight")

    # Set penwidth.  Remove edges that are below threshold

    f_max = 0

    if len(my_weights) > 0:
        f_max = max(my_weights.items(), key=operator.itemgetter(1))[1]

        if not my_args["scale_edge_weight_func"]:
            _scale_edge_weight_func = functools.partial(
                scale_edge_weight,
                f_max=f_max,
                scale=my_args["scale"],
                threshold=my_args["threshold"],
            )
        else:
            _scale_edge_weight_func = functools.partial(
                my_args["scale_edge_weight_func"],
                f_max=f_max,
                scale=my_args["scale"],
                threshold=my_args["threshold"],
            )

        remove_edges = []
        for edge in d_g.edges:
            if not _scale_edge_weight_func(d_g.get_edge_data(*edge)):
                remove_edges.append(edge)

        d_g.remove_edges_from(remove_edges)

    # Remove isolated nodes if desired

    if not my_args["allow_isolated_species"]:
        isolated_nodes = list(nx.isolates(d_g))
        for node in isolated_nodes:
            if (
                node not in my_args["solar_species"]
                and node not in my_args["special_node_attributes"]
            ):
                d_g.remove_node(node)

    # Restore anchors

    for anchor in anchors:
        if anchor not in d_g.nodes:
            d_g.add_node(anchor, style="invis")

    # Get new subset

    sub_graph_2 = nx.subgraph(d_g, subset_nuclides)

    for node in sub_graph_2.nodes:
        sub_graph_2.nodes[node]["pos"] = _get_pos(
            net, node, my_args["state_scaling"]
        )
        sub_graph_2.nodes[node]["label"] = my_args["node_label_func"](node)

    d_g.graph["label"] = my_args["title_func"](f_max)

    _color_edges(sub_graph_2, net, my_args["reaction_color_tuples"])

    return sub_graph_2


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
          `special_node_attributes`_.

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
    ]

    my_args = _get_keywords(my_list, **kwargs)

    result = {}

    subset_nuclides, anchors = _get_subset_and_anchors(
        net.get_nuclides(nuc_xpath=my_args["induced_nuc_xpath"])
    )

    nuc_names = list(net.get_nuclides().keys())

    g_names = net.xml.get_graphviz_names(nuc_names)

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

        result[key] = _create_integrated_current_graph(
            net, value, subset_nuclides, anchors, **my_args
        )

    return result
