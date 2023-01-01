import operator
from math import floor, log10
import networkx as nx
import wnnet.net as wn
import wnnet.zones as wz
import wnnet.flows as wf


def get_solar_species():
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


def make_time_t9_rho_string(time, t9, rho):
    def fexp(f):
        # Add 0.01 for rounding for :.2f mantissa formating
        return int(floor(log10(abs(f)) + 0.01)) if f != 0.0 else 0.0

    def fman(f):
        return f / 10.0 ** fexp(f)

    if time:
        return "<time (s) = {:.2f} x 10<sup>{:d}</sup>, T<sub>9</sub> = {:.2f}, rho (g/cc) = {:.2f} x 10<sup>{:d}</sup> (g/cc)>".format(
            fman(time), fexp(time), t9, fman(rho), fexp(rho)
        )
    else:
        return "<T<sub>9</sub> = {:.2f}, rho (g/cc) = {:.2f} x 10<sup>{:d}</sup> (g/cc)>".format(
            t9, fman(rho), fexp(rho)
        )


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


def create_flow_graph(
    net,
    t9,
    rho,
    mass_fractions,
    time=None,
    induced_nuc_xpath="",
    induced_reac_xpath="",
    reaction_color_tuples=None,
    threshold=0.01,
    node_shape="box",
    scale=10,
    allow_isolated_species=False,
    title_func=None,
    graph_attributes=None,
    edge_attributes=None,
    node_attributes=None,
    solar_species=None,
    solar_node_attributes=None,
    special_node_attributes=None,
):

    result = {}

    nuclides = net.get_nuclides()
    reactions = net.get_reactions()

    f = wf.compute_flows(net, t9, rho, mass_fractions, reac_xpath=induced_reac_xpath)

    # Get the subset of nuclides to view in the graph.  Get anchors.

    val, anchors = _get_subset_and_anchors(net, induced_nuc_xpath)

    DG = nx.MultiDiGraph()

    for nuc in nuclides:
        DG.add_node(nuc, shape=node_shape)

    for r in f:
        tup = f[r]

        if tup[0] > 0:
            for reactant in reactions[r].nuclide_reactants:
                for product in reactions[r].nuclide_products:
                    DG.add_edge(
                        reactant, product, weight=tup[0], reaction=r, arrowsize=0.2
                    )

        if tup[1] > 0:
            for product in reactions[r].nuclide_products:
                for reactant in reactions[r].nuclide_reactants:
                    DG.add_edge(
                        product, reactant, weight=tup[1], reaction=r, arrowsize=0.2
                    )

    # Title

    if title_func:
        DG.graph["label"] = title_func()
    else:
        DG.graph["label"] = make_time_t9_rho_string(time, t9, rho)

    # Apply attributes

    _apply_graph_attributes(DG, graph_attributes)

    _apply_node_attributes(DG, node_attributes)

    _apply_edge_attributes(DG, edge_attributes)

    _apply_solar_node_attributes(DG, solar_species, solar_node_attributes)

    _apply_special_node_attributes(DG, special_node_attributes)

    # Subgraph and maximum flow within subgraph

    S = nx.subgraph(DG, val)

    w = nx.get_edge_attributes(S, "weight")

    f_max = max(w.items(), key=operator.itemgetter(1))[1]

    # Set penwidth.  Remove edges that are below threshold

    remove_edges = []
    for edge in DG.edges:
        penwidth = scale * DG.get_edge_data(*edge)["weight"] / f_max
        if penwidth > threshold:
            DG.get_edge_data(*edge)["penwidth"] = penwidth
        else:
            remove_edges.append(edge)

    DG.remove_edges_from(remove_edges)

    # Remove isolated nodes if desired

    if not allow_isolated_species:
        DG.remove_nodes_from(list(nx.isolates(DG)))

    # Restore anchors

    for anchor in anchors:
        if anchor not in DG.nodes:
            DG.add_node(anchor, style="invis")

    # Get new subset

    S2 = nx.subgraph(DG, val)

    g_names = net.xml.get_graphviz_names(list(S2.nodes.keys()))

    for node in S2.nodes:
        nuc = nuclides[node]
        z = nuc["z"]
        n = nuc["a"] - z
        S2.nodes[node]["pos"] = str(n) + "," + str(z) + "!"
        S2.nodes[node]["label"] = g_names[node]

    _color_edges(S2, net, reaction_color_tuples)

    return S2


def create_zone_flow_graphs(
    net,
    zones,
    induced_nuc_xpath="",
    induced_reac_xpath="",
    reaction_color_tuples=None,
    threshold=0.01,
    node_shape="box",
    scale=10,
    allow_isolated_species=False,
    title_func=None,
    graph_attributes=None,
    edge_attributes=None,
    node_attributes=None,
    solar_species=None,
    solar_node_attributes=None,
    special_node_attributes=None,
):

    result = {}

    # Loop on zones

    for zone in zones:
        s_t9 = "t9"
        s_rho = "rho"
        s_time = "time"
        props = zones[zone]["properties"]
        if s_t9 in props and s_rho in props and s_time in props:
            result[zone] = create_flow_graph(
                net,
                float(props[s_t9]),
                float(props[s_rho]),
                zones[zone]["mass fractions"],
                time=float(props[s_time]),
                induced_nuc_xpath=induced_nuc_xpath,
                induced_reac_xpath=induced_reac_xpath,
                reaction_color_tuples=reaction_color_tuples,
                threshold=threshold,
                node_shape=node_shape,
                scale=scale,
                allow_isolated_species=allow_isolated_species,
                title_func=title_func,
                graph_attributes=graph_attributes,
                edge_attributes=edge_attributes,
                node_attributes=node_attributes,
                solar_species=solar_species,
                solar_node_attributes=solar_node_attributes,
                special_node_attributes=special_node_attributes,
            )

    return result


def create_network_graph(
    net,
    induced_nuc_xpath="",
    induced_reac_xpath="",
    direction="both",
    reaction_color_tuples=None,
    threshold=0.01,
    node_shape="box",
    scale=10,
    allow_isolated_species=False,
    graph_attributes=None,
    edge_attributes=None,
    node_attributes=None,
    solar_species=None,
    solar_node_attributes=None,
    special_node_attributes=None,
):

    assert direction == "forward" or direction == "reverse" or direction == "both"

    result = {}

    nuclides = net.get_nuclides()
    reactions = net.get_reactions(reac_xpath=induced_reac_xpath)

    # Get the subset of nuclides to view in the graph.  Get anchors.

    val, anchors = _get_subset_and_anchors(net, induced_nuc_xpath)

    DG = nx.MultiDiGraph()

    for nuc in nuclides:
        DG.add_node(nuc, shape=node_shape)

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

    g_names = net.xml.get_graphviz_names(list(S.nodes.keys()))

    for node in S.nodes:
        nuc = nuclides[node]
        z = nuc["z"]
        n = nuc["a"] - z
        S.nodes[node]["pos"] = str(n) + "," + str(z) + "!"
        S.nodes[node]["label"] = g_names[node]

    _color_edges(S, net, reaction_color_tuples)

    return S


def create_links_flow_graph(
    net,
    t9,
    rho,
    mass_fractions,
    time=None,
    induced_nuc_xpath="",
    induced_reac_xpath="",
    direction="both",
    reaction_color_tuples=None,
    threshold=0.01,
    node_shape="box",
    scale=10,
    allow_isolated_species=False,
    title_func=None,
    graph_attributes=None,
    edge_attributes=None,
    node_attributes=None,
    solar_species=None,
    solar_node_attributes=None,
    special_node_attributes=None,
):

    result = {}

    nuclides = net.get_nuclides()
    reactions = net.get_reactions()

    f = wf.compute_link_flows(
        net, t9, rho, mass_fractions, reac_xpath=induced_reac_xpath, direction=direction
    )

    # Get the subset of nuclides to view in the graph.  Get anchors.

    val, anchors = _get_subset_and_anchors(net, induced_nuc_xpath)

    DG = nx.MultiDiGraph()

    for nuc in nuclides:
        DG.add_node(nuc, shape=node_shape)

    for r in f:
        tups = f[r]

        for tup in tups:
            if tup[2] > 0:
                DG.add_edge(tup[0], tup[1], weight=tup[2], reaction=r, arrowsize=0.2)

    # Title

    if title_func:
        DG.graph["label"] = title_func()
    else:
        DG.graph["label"] = make_time_t9_rho_string(time, t9, rho)

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


    # Get subset

    S = nx.subgraph(DG, val)

    g_names = net.xml.get_graphviz_names(list(S.nodes.keys()))

    for node in S.nodes:
        nuc = nuclides[node]
        z = nuc["z"]
        n = nuc["a"] - z
        S.nodes[node]["pos"] = str(n) + "," + str(z) + "!"
        S.nodes[node]["label"] = g_names[node]

    _color_edges(S, net, reaction_color_tuples)

    return S
