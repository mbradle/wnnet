"""Module that holds helper routines for the graph module.  This module
is not part of the wnnet API.  This documentation is shown simply for
completeness."""

import operator
from math import floor, log10
import networkx as nx


def get_keywords(kw_list, **kwargs):
    """Helper function to return keyword arguments for graph module."""

    for key in kwargs:
        assert key in kw_list, f"Keyword {key} not allowed in this routine."

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
    def_kwargs["solar_species"] = get_solar_species_list()
    def_kwargs["solar_node_attributes"] = {
        "fillcolor": "yellow",
        "style": "filled",
    }
    def_kwargs["special_node_attributes"] = None
    def_kwargs["special_edge_attributes"] = None

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


def get_solar_species_list():
    """Helper function to return list of solar species."""

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


def color_edges(my_graph, net, color_tuples):
    """Helper function to color arcs in graph by reaction type."""
    color = {}
    for reaction in net.get_reactions():
        color[reaction] = "black"

    if color_tuples:
        for color_tup in color_tuples:
            for reaction in net.get_reactions(reac_xpath=color_tup[0]):
                color[reaction] = color_tup[1]

        for edge in my_graph.edges:
            if not "color" in my_graph.edges[edge]:
                my_graph.edges[edge]["color"] = color[edge[2]]


def fexp(_x):
    """Helper function to compute an exponent of a number."""
    # Add 0.01 for rounding for :.2f mantissa formating
    return int(floor(log10(abs(_x)) + 0.01)) if _x != 0.0 else 0


def fman(_x):
    """Helper function to compute the mantissa of a number."""
    return _x / 10.0 ** fexp(_x)


def get_pos(coll, name, state_scaling):
    """Helper function to define the position of a node in a network graph."""
    _z, _a, state = coll.xml.get_z_a_state_from_nuclide_name(name)
    _n = _a - _z
    if state == "g":
        _z -= state_scaling
    elif state == "m":
        _z += state_scaling
    return (float(_n), float(_z))


def get_subset_and_anchors(nuclides):
    """Helper function to define the scope and anchors of a graph."""
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


def apply_graph_attributes(my_graph, graph_attributes):
    """Helper function to apply attributes to the graph as a whole."""
    if graph_attributes:
        for key, value in graph_attributes.items():
            my_graph.graph[key] = value


def apply_node_attributes(my_graph, node_attributes):
    """Helper function to apply attributes to the nodes of a graph."""
    if node_attributes:
        for key, value in node_attributes.items():
            for node in my_graph.nodes:
                my_graph.nodes[node][key] = value


def apply_edge_attributes(my_graph, edge_attributes):
    """Helper function to apply attributes to the edges of a graph."""
    if edge_attributes:
        for key, value in edge_attributes.items():
            for edge in my_graph.edges:
                my_graph.edges[edge][key] = value


def apply_solar_node_attributes(
    my_graph, solar_species, solar_node_attributes
):
    """Helper function to apply attributes to the nodes representing solar\
       species in a graph."""

    if solar_node_attributes:
        for node in solar_species:
            if node in my_graph.nodes:
                for key, value in solar_node_attributes.items():
                    my_graph.nodes[node][key] = value


def apply_special_node_attributes(my_graph, special_node_attributes):
    """Helper function to apply attributes to special nodes in a graph."""
    if special_node_attributes:
        for node in special_node_attributes:
            for key, value in special_node_attributes[node].items():
                my_graph.nodes[node][key] = value


def apply_special_edge_attributes(my_graph, special_edge_attributes):
    """Helper function to apply attributes to special edges in a graph."""
    if special_edge_attributes:
        for edge in special_edge_attributes:
            if not edge in my_graph.edges:
                my_graph.add_edge(*edge)
            for key, value in special_edge_attributes[edge].items():
                my_graph.edges[edge][key] = value


def create_flow_graph(net, my_flows, subset_nuclides, anchors, **my_args):
    """Helper function to create a flow graph."""

    d_g = nx.MultiDiGraph()

    for nuc in net.get_nuclides():
        d_g.add_node(nuc)

    add_flows_to_graph(net, d_g, my_flows, my_args)

    # Apply attributes

    apply_graph_attributes(d_g, my_args["graph_attributes"])

    apply_node_attributes(d_g, my_args["node_attributes"])

    apply_edge_attributes(d_g, my_args["edge_attributes"])

    apply_solar_node_attributes(
        d_g, my_args["solar_species"], my_args["solar_node_attributes"]
    )

    apply_special_node_attributes(d_g, my_args["special_node_attributes"])

    # Subgraph and widths.

    sub_graph = nx.subgraph(d_g, subset_nuclides)

    f_max = set_widths_and_get_max_weight(d_g, sub_graph, my_args)

    # Apply special edge attributes after scaling.

    apply_special_edge_attributes(d_g, my_args["special_edge_attributes"])

    # Remove isolated nodes if desired

    if not my_args["allow_isolated_species"]:
        remove_isolated_nodes(d_g, my_args)

    # Restore anchors

    for anchor in anchors:
        if anchor not in d_g.nodes:
            d_g.add_node(anchor, style="invis")

    # Get new subset

    sub_graph_2 = nx.subgraph(d_g, subset_nuclides)

    for node in sub_graph_2.nodes:
        sub_graph_2.nodes[node]["pos"] = get_pos(
            net, node, my_args["state_scaling"]
        )
        sub_graph_2.nodes[node]["label"] = my_args["node_label_func"](node)

    # Title

    d_g.graph["label"] = my_args["title_func"](f_max)

    color_edges(sub_graph_2, net, my_args["reaction_color_tuples"])

    return sub_graph_2


def create_integrated_current_graph(
    net, zone, subset_nuclides, anchors, **my_args
):
    """Helper function to create an integrated current graph."""

    # Create graph.

    d_g = nx.MultiDiGraph()

    # Add nodes and arcs.

    for nuc in net.get_nuclides():
        d_g.add_node(nuc)

    add_currents_to_graph(net, zone, d_g)

    # Apply attributes

    apply_graph_attributes(d_g, my_args["graph_attributes"])

    apply_node_attributes(d_g, my_args["node_attributes"])

    apply_edge_attributes(d_g, my_args["edge_attributes"])

    apply_solar_node_attributes(
        d_g, my_args["solar_species"], my_args["solar_node_attributes"]
    )

    apply_special_node_attributes(d_g, my_args["special_node_attributes"])

    # Subgraph and maximum flow within subgraph

    sub_graph = nx.subgraph(d_g, subset_nuclides)

    f_max = set_widths_and_get_max_weight(d_g, sub_graph, my_args)

    # Apply special edge attributes after scaiing.

    apply_special_edge_attributes(d_g, my_args["special_edge_attributes"])

    # Remove isolated nodes if desired

    if not my_args["allow_isolated_species"]:
        remove_isolated_nodes(d_g, my_args)

    # Restore anchors

    for anchor in anchors:
        if anchor not in d_g.nodes:
            d_g.add_node(anchor, style="invis")

    # Get new subset

    sub_graph_2 = nx.subgraph(d_g, subset_nuclides)

    for node in sub_graph_2.nodes:
        sub_graph_2.nodes[node]["pos"] = get_pos(
            net, node, my_args["state_scaling"]
        )
        sub_graph_2.nodes[node]["label"] = my_args["node_label_func"](node)

    d_g.graph["label"] = my_args["title_func"](f_max)

    color_edges(sub_graph_2, net, my_args["reaction_color_tuples"])

    return sub_graph_2


def add_reactions_to_graph(net, d_g, my_args):
    """Helper function to add reaction arcs to graph."""

    reactions = net.get_reactions(reac_xpath=my_args["induced_reac_xpath"])

    for key, value in reactions.items():
        if my_args["direction"] in ("forward", "both"):
            for reactant in value.nuclide_reactants:
                for product in value.nuclide_products:
                    d_g.add_edge(reactant, product, key)

        if not net.is_weak_reaction(key) and (
            my_args["direction"] in ("reverse", "both")
        ):
            for product in value.nuclide_products:
                for reactant in value.nuclide_reactants:
                    d_g.add_edge(product, reactant, key)


def add_flows_to_graph(net, d_g, my_flows, my_args):
    """Helper function to add flow arcs to graph."""

    reactions = net.get_reactions()

    for key, value in my_flows.items():

        if my_args["flow_type"] == "full":

            if value[0] > 0:
                add_flow_edges(
                    d_g,
                    reactions[key].nuclide_reactants,
                    reactions[key].products,
                    value[0],
                    key,
                )

            if value[1] > 0:
                add_flow_edges(
                    d_g,
                    reactions[key].nuclide_products,
                    reactions[key].reactants,
                    value[1],
                    key,
                )

        if my_args["flow_type"] == "net":

            net_flow = value[0] - value[1]

            if net_flow > 0:
                add_flow_edges(
                    d_g,
                    reactions[key].nuclide_reactants,
                    reactions[key].nuclide_products,
                    net_flow,
                    key,
                )
            else:
                add_flow_edges(
                    d_g,
                    reactions[key].nuclide_products,
                    reactions[key].nuclide_reactants,
                    -net_flow,
                    key,
                )


def add_currents_to_graph(net, zone, d_g):
    """Helper function to add current arcs to graph."""

    reactions = net.get_reactions()

    props = zone["properties"]

    f_currents = {}

    for prop in props:
        if isinstance(prop, tuple):
            if prop[0] == "flow current":
                f_currents[prop[1]] = float(props[prop])

    for key, value in f_currents.items():
        if value > 0:
            for reactant in reactions[key].nuclide_reactants:
                for product in reactions[key].nuclide_products:
                    d_g.add_edge(reactant, product, key, weight=value)
        else:
            for product in reactions[key].nuclide_products:
                for reactant in reactions[key].nuclide_reactants:
                    d_g.add_edge(product, reactant, key, weight=-value)


def remove_isolated_nodes(d_g, my_args):
    """Helper function to remove isolated nodes."""

    isolated_nodes = list(nx.isolates(d_g))
    for node in isolated_nodes:
        if node not in my_args["solar_species"] and (
            not my_args["special_node_attributes"]
            or node not in my_args["special_node_attributes"]
        ):
            d_g.remove_node(node)


def set_widths_and_get_max_weight(d_g, sub_graph, my_args):
    """Helper function to set edge widths or remove."""

    my_weights = nx.get_edge_attributes(sub_graph, "weight")

    f_max = 0

    if len(my_weights) > 0:
        f_max = max(my_weights.items(), key=operator.itemgetter(1))[1]

        remove_edges = []
        for edge in d_g.edges:
            if not my_args["scale_edge_weight_func"](
                d_g.get_edge_data(*edge),
                f_max,
                my_args["scale"],
                my_args["threshold"],
            ):
                remove_edges.append(edge)

        d_g.remove_edges_from(remove_edges)

    return f_max


def add_flow_edges(d_g, s_array, t_array, my_weight, key):
    """Helper function to add flow edges."""

    for source in s_array:
        for target in t_array:
            d_g.add_edge(source, target, key, weight=my_weight)
