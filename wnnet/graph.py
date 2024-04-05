import operator
from math import floor, log10
import networkx as nx
import wnnet.net as wn
import wnnet.zones as wz
import wnnet.flows as wf


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


def make_time_t9_rho_flow_string(zone, zone_label, f_max):
    """The default title function for zone flow graphs.

    Args:
        ``zone``:  A `wnutils <https://wnutils.readthedocs.io>`_ zone object.

        ``zone_label`` (:obj:`str` or :obj:`tuple`):  The label of the zone.

        ``f_max`` (:obj:`float`):  The maximum flow in the scope of the graph.
    Returns:
        A :obj:`str`.  The title string.

    """

    def fexp(f):
        # Add 0.01 for rounding for :.2f mantissa formating
        return int(floor(log10(abs(f)) + 0.01)) if f != 0.0 else 0

    def fman(f):
        return f / 10.0 ** fexp(f)

    props = zone["properties"]
    time = float(props["time"])
    t9 = float(props["t9"])
    rho = float(props["rho"])

    return "<time (s) = {:.2f} x 10<sup>{:d}</sup>, T<sub>9</sub> = {:.2f}, rho (g/cc) = {:.2f} x 10<sup>{:d}</sup> (g/cc), max. flow = {:.2f} x 10<sup>{:d}</sup> (s<sup>-1</sup>)>".format(
        fman(time),
        fexp(time),
        t9,
        fman(rho),
        fexp(rho),
        fman(f_max),
        fexp(f_max),
    )


def make_time_t9_rho_current_string(zone, zone_label, f_max):
    """The default title function for zone integrated current graphs.

    Args:
        ``zone``:  A `wnutils <https://wnutils.readthedocs.io>`_ zone object.

        ``zone_label`` (:obj:`str` or :obj:`tuple`):  The label of the zone.

        ``f_max`` (:obj:`float`):  The maximum integrated current in the scope of the graph.
    Returns:
        A :obj:`str`.  The title string.

    """

    def fexp(f):
        # Add 0.01 for rounding for :.2f mantissa formating
        return int(floor(log10(abs(f)) + 0.01)) if f != 0.0 else 0

    def fman(f):
        return f / 10.0 ** fexp(f)

    props = zone["properties"]
    time = float(props["time"])
    t9 = float(props["t9"])
    rho = float(props["rho"])

    return "<time (s) = {:.2f} x 10<sup>{:d}</sup>, T<sub>9</sub> = {:.2f}, rho (g/cc) = {:.2f} x 10<sup>{:d}</sup> (g/cc), max. int. current = {:.2f} x 10<sup>{:d}</sup>>".format(
        fman(time),
        fexp(time),
        t9,
        fman(rho),
        fexp(rho),
        fman(f_max),
        fexp(f_max),
    )


def make_t9_rho_flow_string(f_max, t9, rho):
    """The default title function for flow graphs.

    Args:
        ``f_max`` (:obj:`float`):  The maximum flow in the scope of the graph.

        ``t9`` (:obj:`float`):  The temperature in 10\ :sup:`9` K at which the flows are computed.  This value is bound by the calling function.

        ``rho`` (:obj:`float`):  The density in g/cc at which the flows are computed.  This value is bound by the calling function.

    Returns:
        A :obj:`str`.  The title string.

    """

    def fexp(f):
        # Add 0.01 for rounding for :.2f mantissa formating
        return int(floor(log10(abs(f)) + 0.01)) if f != 0.0 else 0

    def fman(f):
        return f / 10.0 ** fexp(f)

    return "<T<sub>9</sub> = {:.2f}, rho (g/cc) = {:.2f} x 10<sup>{:d}</sup> (g/cc), max. flow = {:.2f} x 10<sup>{:d}</sup> (s<sup>-1</sup>)>".format(
        t9, fman(rho), fexp(rho), fman(f_max), fexp(f_max)
    )


def make_node_label(name, g_names):
    """The default node label function.

    Args:
        ``name`` (:obj:`str`):  The species name.

        ``g_names`` (:obj:`dict`):  A dictionary of graphviz names.  This dictionary is bound to the function by the calling routine.

    Returns:
        A :obj:`str`.  The node label.

    """
    return g_names[name]


def make_zone_node_label(zone, zone_label, name, g_names):
    """The default zone node label function.

    Args:
        ``zone``:  A `wnutils <https://wnutils.readthedocs.io>`_ zone object.

        ``zone_label`` (:obj:`str` or :obj:`tuple`):  The label of the zone.

        ``name`` (:obj:`str`):  The species name.

        ``g_names`` (:obj:`dict`):  A dictionary of graphviz names.  This dictionary is bound to the function by the calling routine.

    Returns:
        A :obj:`str`.  The node label.

    """
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
    r = edge_data["weight"] / f_max
    if r >= threshold:
        edge_data["penwidth"] = scale * r
    else:
        keep_edge = False
    return keep_edge


def _color_edges(G, net, color_tuples):
    color = {}
    for reaction in net.get_reactions():
        color[reaction] = "black"

    if color_tuples:
        for color_tup in color_tuples:
            for reaction in net.get_reactions(reac_xpath=color_tup[0]):
                color[reaction] = color_tup[1]

        for edge in G.edges:
            G.edges[edge]["color"] = color[G.edges[edge]["reaction"]]


def _get_pos(net, name, state_scaling):
    z, a, state = net.xml.get_z_a_state_from_nuclide_name(name)
    n = a - z
    if state == "g":
        z -= state_scaling
    elif state == "m":
        z += state_scaling
    return (n, z)
    # return str(n) + ".," + str(z) + ".!"


def _get_subset_and_anchors(net, induced_nuc_xpath):
    val = []
    dict = {}
    z_dict = {}
    nuclides = net.get_nuclides()
    for sp in net.get_nuclides(nuc_xpath=induced_nuc_xpath):
        val.append(sp)
        dict[(nuclides[sp]["z"], nuclides[sp]["a"])] = sp
        if nuclides[sp]["z"] not in z_dict:
            z_dict[nuclides[sp]["z"]] = []
        else:
            z_dict[nuclides[sp]["z"]].append(nuclides[sp]["a"])

    z_array = []

    for z in z_dict:
        z_array.append(z)
        z_dict[z].sort()

    z_array.sort()

    anchors = []
    anchors.append(dict[(z_array[0], z_dict[z_array[0]][0])])
    anchors.append(dict[(z_array[0], z_dict[z_array[0]][-1])])
    anchors.append(dict[(z_array[-1], z_dict[z_array[-1]][0])])
    anchors.append(dict[(z_array[-1], z_dict[z_array[-1]][-1])])

    return (val, anchors)


def _apply_graph_attributes(G, graph_attributes):
    if graph_attributes:
        for key in graph_attributes:
            G.graph[key] = graph_attributes[key]


def _apply_node_attributes(G, node_attributes):
    if node_attributes:
        for key in node_attributes:
            for node in G.nodes:
                G.nodes[node][key] = node_attributes[key]


def _apply_edge_attributes(G, edge_attributes):
    if edge_attributes:
        for key in edge_attributes:
            for edge in G.edges:
                G.edges[edge][key] = edge_attributes[key]


def _apply_solar_node_attributes(G, solar_species, solar_node_attributes):
    ss_species = solar_species
    ss_node_attributes = solar_node_attributes

    if not solar_species:
        ss_species = get_solar_species()
    if not solar_node_attributes:
        ss_node_attributes = {"fillcolor": "yellow", "style": "filled"}

    for node in ss_species:
        if node in G.nodes:
            for key in ss_node_attributes:
                G.nodes[node][key] = ss_node_attributes[key]


def _apply_special_node_attributes(G, special_node_attributes):
    if special_node_attributes:
        for node in special_node_attributes:
            for key in special_node_attributes[node]:
                G.nodes[node][key] = special_node_attributes[node][key]


def _create_flow_graph(
    net,
    f,
    flow_type,
    subset_nuclides,
    anchors,
    allow_isolated_species,
    reaction_color_tuples,
    threshold,
    scale,
    state_scaling,
    title_func,
    node_label_func,
    scale_edge_weight_func,
    graph_attributes,
    node_attributes,
    edge_attributes,
    solar_species,
    solar_node_attributes,
    special_node_attributes,
):

    nuclides = net.get_nuclides()
    reactions = net.get_reactions()

    # Solar species

    _solar_species = solar_species
    if not solar_species:
        _solar_species = get_solar_species()

    # Solar species

    _solar_species = solar_species
    if not solar_species:
        _solar_species = get_solar_species()

    DG = nx.MultiDiGraph()

    for nuc in nuclides:
        DG.add_node(nuc, shape="box", fontsize=16)

    for r in f:
        tup = f[r]

        if flow_type == "full":
            if tup[0] > 0:
                for reactant in reactions[r].nuclide_reactants:
                    for product in reactions[r].nuclide_products:
                        DG.add_edge(
                            reactant,
                            product,
                            weight=tup[0],
                            reaction=r,
                            arrowsize=0.2,
                        )

            if tup[1] > 0:
                for product in reactions[r].nuclide_products:
                    for reactant in reactions[r].nuclide_reactants:
                        DG.add_edge(
                            product,
                            reactant,
                            weight=tup[1],
                            reaction=r,
                            arrowsize=0.2,
                        )

        elif flow_type == "net":
            net_flow = tup[0] - tup[1]

            if net_flow > 0:
                for reactant in reactions[r].nuclide_reactants:
                    for product in reactions[r].nuclide_products:
                        DG.add_edge(
                            reactant,
                            product,
                            weight=net_flow,
                            reaction=r,
                            arrowsize=0.2,
                        )

            if net_flow < 0:
                for product in reactions[r].nuclide_products:
                    for reactant in reactions[r].nuclide_reactants:
                        DG.add_edge(
                            product,
                            reactant,
                            weight=-net_flow,
                            reaction=r,
                            arrowsize=0.2,
                        )

    # Apply attributes

    _apply_graph_attributes(DG, graph_attributes)

    _apply_node_attributes(DG, node_attributes)

    _apply_edge_attributes(DG, edge_attributes)

    _apply_solar_node_attributes(DG, _solar_species, solar_node_attributes)

    _apply_special_node_attributes(DG, special_node_attributes)

    # Subgraph and maximum flow within subgraph

    S = nx.subgraph(DG, subset_nuclides)

    w = nx.get_edge_attributes(S, "weight")

    # Set penwidth.  Remove edges that are below threshold

    f_max = 0

    if len(w) > 0:
        f_max = max(w.items(), key=operator.itemgetter(1))[1]

        if not scale_edge_weight_func:
            _scale_edge_weight_func = lambda edge_data: scale_edge_weight(
                edge_data, f_max, scale, threshold
            )
        else:
            _scale_edge_weight_func = lambda edge_data: scale_edge_weight_func(
                edge_data, f_max, scale, threshold
            )

        remove_edges = []
        for edge in DG.edges:
            if not _scale_edge_weight_func(DG.get_edge_data(*edge)):
                remove_edges.append(edge)

        DG.remove_edges_from(remove_edges)

    # Remove isolated nodes if desired

    if not allow_isolated_species:
        isolated_nodes = list(nx.isolates(DG))
        for node in isolated_nodes:
            if node not in _solar_species:
                DG.remove_node(node)

    # Restore anchors

    for anchor in anchors:
        if anchor not in DG.nodes:
            DG.add_node(anchor, style="invis")

    # Get new subset

    S2 = nx.subgraph(DG, subset_nuclides)

    for node in S2.nodes:
        S2.nodes[node]["pos"] = _get_pos(net, node, state_scaling)
        S2.nodes[node]["label"] = node_label_func(node)

    # Title

    DG.graph["label"] = title_func(f_max)

    _color_edges(S2, net, reaction_color_tuples)

    return S2


def create_flow_graph(
    net,
    t9,
    rho,
    mass_fractions,
    flow_type="net",
    induced_nuc_xpath="",
    induced_reac_xpath="",
    user_funcs="",
    reaction_color_tuples=None,
    threshold=0.01,
    scale=10,
    state_scaling=0.325,
    allow_isolated_species=False,
    title_func=None,
    node_label_func=None,
    scale_edge_weight_func=None,
    graph_attributes=None,
    edge_attributes=None,
    node_attributes=None,
    solar_species=None,
    solar_node_attributes=None,
    special_node_attributes=None,
):
    """A routine to create a flow graph for a given set of mass fractions at the input temperature and density.

    Args:
        ``net``: A wnnet network. 

        ``t9`` (:obj:`float`):  The temperature in 10\ :sup:`9` K at which to compute the flows. 

        ``rho`` (:obj:`float`):  The density in g/cc at which to compute the flows.
        
        ``mass_fractions`` (:obj:`float`): A `wnutils <https://wnutils.readthedocs.io>`_ dictionary of mass fractions.

        ``flow_type`` (:obj:`str`, optional): A string giving the flow type to be presented.  The possible values are `net`, which shows the forward minus the reverse flow (or the opposite if the reverse flow is larger), and `full`, which shows both the foward and reverse flows.

        ``induced_nuc_xpath`` (:obj:`str`, optional): An XPath expression to select the subset of nuclides in the graph.  The default is all species in the network.

        ``induced_reac_xpath`` (:obj:`str`, optional): An XPath expression to select the subset of reactions in the graph.  The default is all reactions in the network.

        ``reaction_color_tuples`` (:obj:`tuple`, optional): A tuple to select arc colors for reaction types.  The first member of the tuple is an XPath expression to select the reaction type while the second member is a string giving the color for that reaction type.  The default is that all arcs are black.

        ``user_funcs`` (:obj:`dict`, optional): A dictionary of user-defined
        functions associated with a user_rate key.
        The prototype for each
        user rate function should be (*reaction*, *t9*), where
        *t9* is the temperature in billions of Kelvin and *reaction*
        is a `wnutils <https://wnutils.readthedocs.io>`_ reaction
        instance.  Other data can be bound to the function.

        ``threshold`` (:obj:`float`, optional):  The minimum flow (relative to the maximum flow) to be shown on the graph
        
        ``scale`` (:obj:`float`, optional):  Scaling factor for the maximum weight arc.
        
        ``state_scaling`` (:obj:`float`, optional):  Scaling factor for isomeric states.
        
        ``allow_isolated_species`` (:obj:`bool`, optional):  Boolean to choose whether to allow isolated species (ones without incoming or outgoing arcs) in the graph.

        ``title_func`` (optional): A `function <https://docs.python.org/3/library/stdtypes.html#functions>`_ that applies the title to the graph.  The function must take one :obj:`float` argument giving the maximum flow. Other data can be bound to the function.  The function must return a :obj:`str` giving the title.  The default is :meth:`wnnet.graph.make_t9_rho_flow_string`.
        
        ``node_label_func`` (optional): A `function \
            <https://docs.python.org/3/library/stdtypes.html#functions>`_ \
            that applies a label to each node in the graph.  The function \
            must take as argument a species name.  Other data can be bound to \
            the function.  The function must return a :obj:`str` \
            giving the label.  \
            The default is :meth:`wnnet.graph.make_node_label`.
        
        ``scale_edge_weight_func`` (optional): A `function \
            <https://docs.python.org/3/library/stdtypes.html#functions>`_ \
            that applies scales each edge weight in the graph.  The function \
            must take as four arguments: a dictionary of edge data, the \
            maximum edge weight in the scope of the graph, a scale factor \
            by which to scale the weight (input as *scale* to this routine), \
            and a threshold for not including the edge in the graph \
            (input as *threshold* to this routine). \
            Other data can be bound to \
            the function.  The function must modify the weight in the \
            edge data and return a :obj:`bool` indicating whether to include \
            the edge in the graph (True) or not (False).\
            The default is :meth:`wnnet.graph.scale_edge_weight`.
        
        ``graph_attributes`` (:obj:`dict`, optional):  A dictionary of graphviz attributes for the graph.

        ``edge_attributes`` (:obj:`dict`, optional):  A dictionary of grapvhiz attributes for the edges in the graph.

        ``node_attributes`` (:obj:`dict`, optional):  A dictionary of graphviz attributes for the nodes in the graph.

        ``solar_species`` (:obj:`list`, optional):  A list of species to be considered as the naturally occurring species.  The default is the list returned from :meth:`wnnet.graph.get_solar_species`.

        ``solar_node_attributes`` (:obj:`dict`, optional):  A dictionary of graphviz attributes to be applied to the solar species in the graph.

        ``special_node_attributes`` (:obj:`dict`, optional):  A dictionary of graphviz attributes to be applied to the special nodes in the graph.  The dictionary has as keys the names of the special nodes and as values a dictionary of graphviz properties to be applied to the given special node.

    Returns:
        A `networkx multidigraph <https://networkx.org/documentation/stable/reference/classes/multidigraph.html>`_ showing the flows.
     

    """
    assert flow_type == "net" or flow_type == "full"

    f = wf.compute_flows(
        net,
        t9,
        rho,
        mass_fractions,
        reac_xpath=induced_reac_xpath,
        user_funcs=user_funcs,
    )

    # Get the subset of nuclides to view in the graph.  Get anchors.

    subset_nuclides, anchors = _get_subset_and_anchors(net, induced_nuc_xpath)

    # Title

    if not title_func:
        _title_func = lambda f_max: make_t9_rho_flow_string(f_max, t9, rho)
    else:
        _title_func = title_func

    # Node label

    if not node_label_func:
        g_names = net.xml.get_graphviz_names(subset_nuclides)
        _node_label_func = lambda name: make_node_label(name, g_names)
    else:
        _node_label_func = node_label_func

    return _create_flow_graph(
        net,
        f,
        flow_type,
        subset_nuclides,
        anchors,
        allow_isolated_species,
        reaction_color_tuples,
        threshold,
        scale,
        state_scaling,
        _title_func,
        _node_label_func,
        scale_edge_weight_func,
        graph_attributes,
        node_attributes,
        edge_attributes,
        solar_species,
        solar_node_attributes,
        special_node_attributes,
    )


def create_zone_flow_graphs(
    net,
    zones,
    flow_type="net",
    induced_nuc_xpath="",
    induced_reac_xpath="",
    reaction_color_tuples=None,
    user_funcs="",
    threshold=0.01,
    scale=10,
    state_scaling=0.325,
    allow_isolated_species=False,
    title_func=None,
    zone_node_label_func=None,
    scale_edge_weight_func=None,
    graph_attributes=None,
    edge_attributes=None,
    node_attributes=None,
    solar_species=None,
    solar_node_attributes=None,
    special_node_attributes=None,
):
    """A routine to create flow graphs for a set of zones.

    Args:
        ``net``: A wnnet network. 

        ``zones``: A `wnutils <https://wnutils.readthedocs.io>`_ dictionary of zones.

        ``flow_type`` (:obj:`str`, optional): A string giving the flow type to be presented.  The possible values are `net`, which shows the forward minus the reverse flow (or the opposite if the reverse flow is larger), and `full`, which shows both the foward and reverse flows.

        ``induced_nuc_xpath`` (:obj:`str`, optional): An XPath expression to select the subset of nuclides in the graph.  The default is all species in the network.

        ``induced_reac_xpath`` (:obj:`str`, optional): An XPath expression to select the subset of reactions in the graph.  The default is all reactions in the network.

        ``user_funcs`` (:obj:`dict`, optional): A dictionary of user-defined
        functions associated with a user_rate key.
        The prototype for each
        user rate function should be (*reaction*, *t9*, *zone*), where
        *t9* is the temperature in billions of Kelvin and *reaction* and
        *zone* are `wnutils <https://wnutils.readthedocs.io>`_ reaction and
        zone instances.  Other data can be bound to the function.

        ``reaction_color_tuples`` (:obj:`tuple`, optional): A tuple to select arc colors for reaction types.  The first member of the tuple is an XPath expression to select the reaction type while the second member is a string giving the color for that reaction type.  The default is that all arcs are black.

        ``threshold`` (:obj:`float`, optional):  The minimum flow (relative to the maximum flow) to be shown on the graph
        
        ``scale`` (:obj:`float`, optional):  Scaling factor for the maximum weight arc.
        
        ``state_scaling`` (:obj:`float`, optional):  Scaling factor for isomeric states.
        
        ``allow_isolated_species`` (:obj:`bool`, optional):  Boolean to choose whether to allow isolated species (ones without incoming or outgoing arcs) in the graph.

        ``title_func`` (optional): A `function \
             <https://docs.python.org/3/library/stdtypes.html#functions>`_ \
             that applies the title to the graph.  The function must take \
             three arguments.  The first is the zone object corresponding to\
             the graph while the second is the zone label and the third is \
             a :obj:`float` giving the maximum flow. Other data can \
             be bound to the function.  The function must return a :obj:`str` \
             giving the title.  \
             The default is :meth:`wnnet.graph.make_time_t9_rho_flow_string.`.
        
        ``zone_node_label_func`` (optional): A `function \
            <https://docs.python.org/3/library/stdtypes.html#functions>`_ \
            that applies a label to each node in the graph.  The function \
            must take as arguments a zone, the zone label, and a species name. \
            Other data can be bound to \
            the function.  The function must return a :obj:`str` \
            giving the label.  The default is \
            :meth:`wnnet.graph.make_zone_node_label`.
        
        ``scale_edge_weight_func`` (optional): A `function \
            <https://docs.python.org/3/library/stdtypes.html#functions>`_ \
            that applies scales each edge weight in the graph.  The function \
            must take as four arguments: a dictionary of edge data, the \
            maximum edge weight in the scope of the graph, a scale factor \
            by which to scale the weight (input as *scale* to this routine), \
            and a threshold for not including the edge in the graph \
            (input as *threshold* to this routine). \
            Other data can be bound to \
            the function.  The function must modify the weight in the \
            edge data and return a :obj:`bool` indicating whether to include \
            the edge in the graph (True) or not (False).\
            The default is :meth:`wnnet.graph.scale_edge_weight`.
        
        ``graph_attributes`` (:obj:`dict`, optional):  A dictionary of graphviz attributes for the graph.

        ``edge_attributes`` (:obj:`dict`, optional):  A dictionary of grapvhiz attributes for the edges in the graph.

        ``node_attributes`` (:obj:`dict`, optional):  A dictionary of graphviz attributes for the nodes in the graph.

        ``solar_species`` (:obj:`list`, optional):  A list of species to be considered as the naturally occurring species.  The default is the list returned from :meth:`wnnet.graph.get_solar_species`.

        ``solar_node_attributes`` (:obj:`dict`, optional):  A dictionary of graphviz attributes to be applied to the solar species in the graph.

        ``special_node_attributes`` (:obj:`dict`, optional):  A dictionary of graphviz attributes to be applied to the special nodes in the graph.  The dictionary has as keys the names of the special nodes and as values a dictionary of graphviz properties to be applied to the given special node.

    Returns:
        A :obj:`dict` of `networkx multidigraphs <https://networkx.org/documentation/stable/reference/classes/multidigraph.html>`_ showing the flows.  The keys are the zone labels.
     

    """

    result = {}

    f = wf.compute_flows_for_zones(
        net, zones, reac_xpath=induced_reac_xpath, user_funcs=user_funcs
    )

    subset_nuclides, anchors = _get_subset_and_anchors(net, induced_nuc_xpath)

    # Loop on zones

    for zone in f:

        # Title

        if not title_func:
            _title_func = lambda f_max: make_time_t9_rho_flow_string(
                zones[zone], zone, f_max
            )
        else:
            _title_func = lambda f_max: title_func(zones[zone], zone, f_max)

        # Node label

        g_names = net.xml.get_graphviz_names(subset_nuclides)
        if not zone_node_label_func:
            _zone_node_label_func = lambda name: make_zone_node_label(
                zones[zone], zone, name, g_names
            )
        else:
            _zone_node_label_func = lambda name: zone_node_label_func(
                zones[zone], zone, name
            )

        # Create graph

        result[zone] = _create_flow_graph(
            net,
            f[zone],
            flow_type,
            subset_nuclides,
            anchors,
            allow_isolated_species,
            reaction_color_tuples,
            threshold,
            scale,
            state_scaling,
            _title_func,
            _zone_node_label_func,
            scale_edge_weight_func,
            graph_attributes,
            node_attributes,
            edge_attributes,
            solar_species,
            solar_node_attributes,
            special_node_attributes,
        )

    return result


def create_network_graph(
    net,
    induced_nuc_xpath="",
    induced_reac_xpath="",
    direction="both",
    reaction_color_tuples=None,
    threshold=0.01,
    scale=10,
    state_scaling=0.325,
    allow_isolated_species=False,
    node_label_func=None,
    graph_attributes=None,
    edge_attributes=None,
    node_attributes=None,
    solar_species=None,
    solar_node_attributes=None,
    special_node_attributes=None,
):
    """A routine to create a network graph showing species and reactions among them.

    Args:
        ``net``: A wnnet network. 

        ``induced_nuc_xpath`` (:obj:`str`, optional): An XPath expression to select the subset of nuclides in the graph.  The default is all species in the network.

        ``induced_reac_xpath`` (:obj:`str`, optional): An XPath expression to select the subset of reactions in the graph.  The default is all reactions in the network.

        ``direction`` (:obj:`str`, optional): A string indicting the reaction directions to show.  Allowed values are `forward`, `reverse`, and `both` (the default, which shows both `forward` and `reverse`).

        ``reaction_color_tuples`` (:obj:`tuple`, optional): A tuple to select arc colors for reaction types.  The first member of the tuple is an XPath expression to select the reaction type while the second member is a string giving the color for that reaction type.  The default is that all arcs are black.

        ``threshold`` (:obj:`float`, optional):  The minimum flow (relative to the maximum flow) to be shown on the graph
        
        ``scale`` (:obj:`float`, optional):  Scaling factor for the maximum weight arc.
        
        ``state_scaling`` (:obj:`float`, optional):  Scaling factor for isomeric states.
        
        ``allow_isolated_species`` (:obj:`bool`, optional):  Boolean to choose whether to allow isolated species (ones without incoming or outgoing arcs) in the graph.

        ``node_label_func`` (optional): A `function \
            <https://docs.python.org/3/library/stdtypes.html#functions>`_ \
            that applies a label to each node in the graph.  The function \
            must take as argument a species name.  Other data can be bound to \
            the function.  The function must return a :obj:`str` \
            giving the label.  The default is \
            :meth:`wnnet.graph.make_node_label`.
        
        ``graph_attributes`` (:obj:`dict`, optional):  A dictionary of graphviz attributes for the graph.

        ``edge_attributes`` (:obj:`dict`, optional):  A dictionary of grapvhiz attributes for the edges in the graph.

        ``node_attributes`` (:obj:`dict`, optional):  A dictionary of graphviz attributes for the nodes in the graph.

        ``solar_species`` (:obj:`list`, optional):  A list of species to be considered as the naturally occurring species.  The default is the list returned from :meth:`wnnet.graph.get_solar_species`.

        ``solar_node_attributes`` (:obj:`dict`, optional):  A dictionary of graphviz attributes to be applied to the solar species in the graph.

        ``special_node_attributes`` (:obj:`dict`, optional):  A dictionary of graphviz attributes to be applied to the special nodes in the graph.  The dictionary has as keys the names of the special nodes and as values a dictionary of graphviz properties to be applied to the given special node.

    Returns:
        A `networkx multidigraph <https://networkx.org/documentation/stable/reference/classes/multidigraph.html>`_ showing the network; that is, the nuclides and reactions among them.
     
    """

    assert (
        direction == "forward" or direction == "reverse" or direction == "both"
    )

    result = {}

    nuclides = net.get_nuclides()
    reactions = net.get_reactions(reac_xpath=induced_reac_xpath)

    # Get the subset of nuclides to view in the graph.  Get anchors.

    val, anchors = _get_subset_and_anchors(net, induced_nuc_xpath)

    DG = nx.MultiDiGraph()

    for nuc in nuclides:
        DG.add_node(nuc, shape="box", fontsize=16)

    for r in reactions:
        if direction == "forward" or direction == "both":
            for reactant in reactions[r].nuclide_reactants:
                for product in reactions[r].nuclide_products:
                    DG.add_edge(reactant, product, reaction=r, arrowsize=0.2)

        if not net.is_weak_reaction(r) and (
            direction == "reverse" or direction == "both"
        ):
            for product in reactions[r].nuclide_products:
                for reactant in reactions[r].nuclide_reactants:
                    DG.add_edge(product, reactant, reaction=r, arrowsize=0.2)

    # Apply attributes

    _apply_graph_attributes(DG, graph_attributes)

    _apply_node_attributes(DG, node_attributes)

    _apply_edge_attributes(DG, edge_attributes)

    _apply_solar_node_attributes(DG, solar_species, solar_node_attributes)

    _apply_special_node_attributes(DG, special_node_attributes)

    # Remove isolated nodes if desired

    if not allow_isolated_species:
        DG.remove_nodes_from(list(nx.isolates(DG)))

    # Restore anchors

    for anchor in anchors:
        if anchor not in DG.nodes:
            DG.add_node(anchor, style="invis")

    # Subgraph

    S = nx.subgraph(DG, val)

    # Node label

    if node_label_func:
        _node_label_func = node_label_func
    else:
        g_names = net.xml.get_graphviz_names(list(S.nodes.keys()))
        _node_label_func = lambda name: make_node_label(name, g_names)

    for node in S.nodes:
        S.nodes[node]["pos"] = _get_pos(net, node, state_scaling)
        S.nodes[node]["label"] = _node_label_func(node)

    _color_edges(S, net, reaction_color_tuples)

    return S


def _create_integrated_current_graph(
    net,
    zone,
    subset_nuclides,
    induced_nuc_xpath,
    induced_reac_xpath,
    reaction_color_tuples,
    threshold,
    scale,
    state_scaling,
    allow_isolated_species,
    title_func,
    zone_node_label_func,
    scale_edge_weight_func,
    graph_attributes,
    edge_attributes,
    node_attributes,
    solar_species,
    solar_node_attributes,
    special_node_attributes,
):
    nuclides = net.get_nuclides()
    reactions = net.get_reactions()

    # Solar species

    _solar_species = solar_species
    if not solar_species:
        _solar_species = get_solar_species()

    # Solar species

    _solar_species = solar_species
    if not solar_species:
        _solar_species = get_solar_species()

    subset_nuclides, anchors = _get_subset_and_anchors(net, induced_nuc_xpath)

    props = zone["properties"]

    DG = nx.MultiDiGraph()

    for nuc in nuclides:
        DG.add_node(nuc, shape="box", fontsize=16)

    f = {}

    for prop in props:
        if isinstance(prop, tuple):
            if prop[0] == "flow current":
                f[prop[1]] = float(props[prop])

    for r in f:

        if f[r] > 0:
            for reactant in reactions[r].nuclide_reactants:
                for product in reactions[r].nuclide_products:
                    DG.add_edge(
                        reactant,
                        product,
                        weight=f[r],
                        reaction=r,
                        arrowsize=0.2,
                    )

        if f[r] < 0:
            for product in reactions[r].nuclide_products:
                for reactant in reactions[r].nuclide_reactants:
                    DG.add_edge(
                        product,
                        reactant,
                        weight=-f[r],
                        reaction=r,
                        arrowsize=0.2,
                    )

    # Apply attributes

    _apply_graph_attributes(DG, graph_attributes)

    _apply_node_attributes(DG, node_attributes)

    _apply_edge_attributes(DG, edge_attributes)

    _apply_solar_node_attributes(DG, _solar_species, solar_node_attributes)

    _apply_special_node_attributes(DG, special_node_attributes)

    # Subgraph and maximum flow within subgraph

    S = nx.subgraph(DG, subset_nuclides)

    w = nx.get_edge_attributes(S, "weight")

    # Set penwidth.  Remove edges that are below threshold

    f_max = 0

    if len(w) > 0:
        f_max = max(w.items(), key=operator.itemgetter(1))[1]

        if not scale_edge_weight_func:
            _scale_edge_weight_func = lambda edge_data: scale_edge_weight(
                edge_data, f_max, scale, threshold
            )
        else:
            _scale_edge_weight_func = lambda edge_data: scale_edge_weight_func(
                edge_data, f_max, scale, threshold
            )

        remove_edges = []
        for edge in DG.edges:
            if not _scale_edge_weight_func(DG.get_edge_data(*edge)):
                remove_edges.append(edge)

        DG.remove_edges_from(remove_edges)

    # Remove isolated nodes if desired

    if not allow_isolated_species:
        isolated_nodes = list(nx.isolates(DG))
        for node in isolated_nodes:
            if node not in _solar_species:
                DG.remove_node(node)

    # Restore anchors

    for anchor in anchors:
        if anchor not in DG.nodes:
            DG.add_node(anchor, style="invis")

    # Get new subset

    S2 = nx.subgraph(DG, subset_nuclides)

    for node in S2.nodes:
        S2.nodes[node]["pos"] = _get_pos(net, node, state_scaling)
        S2.nodes[node]["label"] = zone_node_label_func(node)

    DG.graph["label"] = title_func(f_max)

    _color_edges(S2, net, reaction_color_tuples)

    return S2


def create_zone_integrated_current_graphs(
    net,
    zones,
    flow_type="net",
    induced_nuc_xpath="",
    induced_reac_xpath="",
    reaction_color_tuples=None,
    threshold=0.01,
    scale=10,
    state_scaling=0.325,
    allow_isolated_species=False,
    title_func=None,
    zone_node_label_func=None,
    scale_edge_weight_func=None,
    graph_attributes=None,
    edge_attributes=None,
    node_attributes=None,
    solar_species=None,
    solar_node_attributes=None,
    special_node_attributes=None,
):
    """A routine to create flow graphs for a set of zones.

    Args:
        ``net``: A wnnet network. 

        ``zones``: A `wnutils <https://wnutils.readthedocs.io>`_ dictionary of zones.

        ``flow_type`` (:obj:`str`, optional): A string giving the flow type to be presented.  The possible values are `net`, which shows the forward minus the reverse flow (or the opposite if the reverse flow is larger), and `full`, which shows both the foward and reverse flows.

        ``induced_nuc_xpath`` (:obj:`str`, optional): An XPath expression to select the subset of nuclides in the graph.  The default is all species in the network.

        ``induced_reac_xpath`` (:obj:`str`, optional): An XPath expression to select the subset of reactions in the graph.  The default is all reactions in the network.

        ``reaction_color_tuples`` (:obj:`tuple`, optional): A tuple to select arc colors for reaction types.  The first member of the tuple is an XPath expression to select the reaction type while the second member is a string giving the color for that reaction type.  The default is that all arcs are black.

        ``threshold`` (:obj:`float`, optional):  The minimum flow (relative to the maximum flow) to be shown on the graph
        
        ``scale`` (:obj:`float`, optional):  Scaling factor for the maximum weight arc.
        
        ``state_scaling`` (:obj:`float`, optional):  Scaling factor for isomeric states.
        
        ``allow_isolated_species`` (:obj:`bool`, optional):  Boolean to choose whether to allow isolated species (ones without incoming or outgoing arcs) in the graph.

        ``title_func`` (optional): A `function \
             <https://docs.python.org/3/library/stdtypes.html#functions>`_ \
             that applies the title to the graph.  The function must take \
             three arguments.  The first is the zone object corresponding to \
             the graph while the second is the zone label (a string or tuple) \
             and the third is a :obj:`float` giving the maximum \
             integrated current.  Other data can \
             be bound to the function.  The function must return a :obj:`str` \
             giving the title.  The default is \
             :meth:`wnnet.graph.make_time_t9_rho_current_string`.
        
        ``zone_node_label_func`` (optional): A `function \
            <https://docs.python.org/3/library/stdtypes.html#functions>`_ \
            that applies a label to each node in the graph.  The function \
            must take three arguments: a zone, its label, and a species name. \
            Other data can be bound to \
            the function.  The function must return a :obj:`str` \
            giving the label.  The default is \
            :meth:`wnnet.graph.make_zone_node_label`.
        
        ``scale_edge_weight_func`` (optional): A `function \
            <https://docs.python.org/3/library/stdtypes.html#functions>`_ \
            that applies scales each edge weight in the graph.  The function \
            must take as four arguments: a dictionary of edge data, the \
            maximum edge weight in the scope of the graph, a scale factor \
            by which to scale the weight (input as *scale* to this routine), \
            and a threshold for not including the edge in the graph \
            (input as *threshold* to this routine). \
            Other data can be bound to \
            the function.  The function must modify the weight in the \
            edge data and return a :obj:`bool` indicating whether to include \
            the edge in the graph (True) or not (False).\
            The default is :meth:`wnnet.graph.scale_edge_weight`.
        
        ``graph_attributes`` (:obj:`dict`, optional):  A dictionary of graphviz attributes for the graph.

        ``edge_attributes`` (:obj:`dict`, optional):  A dictionary of grapvhiz attributes for the edges in the graph.

        ``node_attributes`` (:obj:`dict`, optional):  A dictionary of graphviz attributes for the nodes in the graph.

        ``solar_species`` (:obj:`list`, optional):  A list of species to be considered as the naturally occurring species.  The default is the Solar System's actual naturally occurring species.

        ``solar_node_attributes`` (:obj:`dict`, optional):  A dictionary of graphviz attributes to be applied to the solar species in the graph.

        ``special_node_attributes`` (:obj:`dict`, optional):  A dictionary of graphviz attributes to be applied to the special nodes in the graph.  The dictionary has as keys the names of the special nodes and as values a dictionary of graphviz properties to be applied to the given special node.

    Returns:
        A :obj:`dict` of `networkx multidigraphs <https://networkx.org/documentation/stable/reference/classes/multidigraph.html>`_ showing the integrated currents.  The keys are the zone labels.

    """

    result = {}

    subset_nuclides, anchors = _get_subset_and_anchors(net, induced_nuc_xpath)

    for zone in zones:

        # Title

        if not title_func:
            _title_func = lambda f_max: make_time_t9_rho_current_string(
                zones[zone], zone, f_max
            )
        else:
            _title_func = lambda f_max: title_func(zones[zone], zone, f_max)

        # Node label

        g_names = net.xml.get_graphviz_names(subset_nuclides)
        if not zone_node_label_func:
            _zone_node_label_func = lambda name: make_zone_node_label(
                zones[zone], zone, name, g_names
            )
        else:
            _zone_node_label_func = lambda name: zone_node_label_func(
                zones[zone], zone, name
            )

        result[zone] = _create_integrated_current_graph(
            net,
            zones[zone],
            subset_nuclides,
            induced_nuc_xpath,
            induced_reac_xpath,
            reaction_color_tuples,
            threshold,
            scale,
            state_scaling,
            allow_isolated_species,
            _title_func,
            _zone_node_label_func,
            scale_edge_weight_func,
            graph_attributes,
            edge_attributes,
            node_attributes,
            solar_species,
            solar_node_attributes,
            special_node_attributes,
        )

    return result
