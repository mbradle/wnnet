import requests, io
import wnnet as wn


def get_net():
    return wn.net.Net(
        io.BytesIO(requests.get("https://osf.io/4gmyr/download").content)
    )


def get_nuc():
    return wn.nuc.Nuc(
        io.BytesIO(requests.get("https://osf.io/4gmyr/download").content)
    )


def get_zones():
    return wn.zones.Zones(
        io.BytesIO(requests.get("https://osf.io/4gmyr/download").content)
    )


def test_nuc():
    net = get_net()
    nuc = net.get_nuclides()
    assert len(nuc) > 0
    assert net.compute_nuclear_partition_function("ca41", 6.7) > 0
    assert net.compute_quantum_abundance("o16", 0.2, 10) > 0
    assert net.compute_binding_energy("ca40") > 0


def test_reac():
    net = get_net()
    reac = net.get_reactions()
    assert net.is_weak_reaction("nd135 -> pr135 + positron + neutrino_e")
    assert not net.is_weak_reaction("h1 + cd101 -> he4 + ag98")
    reac_ng = net.get_reactions(
        reac_xpath="[reactant = 'n' and product = 'gamma']"
    )
    assert len(reac) > len(reac_ng) > 0
    for r in reac_ng.values():
        assert r.compute_rate(1.0) > 0


def test_net():
    net = get_net()
    q = net.compute_q_values()
    assert q["n + v46 -> v47 + gamma"] > 0
    rates = net.compute_rates(2)
    assert (
        rates["n + v46 -> v47 + gamma"][0] > 0
        and rates["n + v46 -> v47 + gamma"][1] > 0
    )


def test_flows():
    net = get_net()
    zn = get_zones()
    zones = zn.get_zones(zone_xpath="[position() >= 20 and position() <= 25]")

    t_9 = 1
    rho = 1.0e4
    mass_fractions = {
        ("n", 0, 1): 0.2,
        ("h1", 1, 1): 0.3,
        ("he4", 2, 4): 0.1,
        ("c12", 6, 12): 0.1,
        ("o16", 8, 16): 0.1,
        ("fe56", 26, 56): 0.2,
    }

    f = wn.flows.compute_flows(net, t_9, rho, mass_fractions)
    for val in f.values():
        assert len(val) > 0

    f = wn.flows.compute_flows_for_zones(net, zones)
    for val in f.values():
        assert len(val) > 0

    f = wn.flows.compute_link_flows_for_zones(net, zones, include_dt=True)
    for val in f.values():
        assert len(val) > 0


def test_graph():
    nuc = get_nuc()
    net = get_net()
    zn = get_zones()
    zones = zn.get_zones(zone_xpath="[position() >= 20 and position() <= 25]")

    my_induced_nuc_xpath = "[z <= 30]"
    my_induced_reac_xpath = "[(reactant = 'h1' and product = 'gamma') or\
        (product = 'electron') or (reactant = 'electron') or \
        (product = 'positron')]"
    my_graph_attributes = {
        "label": "My Graph",
        "labelloc": "t",
        "fontsize": 30,
        "fontcolor": "blue",
    }
    my_edge_attributes = {"penwidth": 2}
    my_node_attributes = {"style": "filled", "fillcolor": "bisque"}
    my_special_node_attributes = {
        "fe56": {"fillcolor": "green", "shape": "oval", "style": "filled"}
    }
    my_special_edge_attributes = {
        ("fe56", "fe57", "n + fe56 -> fe57 + gamma"): {
            "color": "blue",
            "penwidth": 2,
        },
        ("co51", "fe51", "n + co51 -> h1 + fe51"): {
            "color": "orange",
            "penwidth": 2,
        },
        ("co51", "fe51", "o51 -> fe51 + positron + neutrino_e"): {
            "color": "violet",
            "penwidth": 2,
        },
    }

    my_color_tuples = [
        ("[product = 'electron']", "blue"),
        ("[(reactant = 'electron') or (product = 'positron')]", "red"),
    ]

    assert wn.graph.create_network_graph(
        net,
        induced_nuc_xpath=my_induced_nuc_xpath,
        induced_reac_xpath=my_induced_reac_xpath,
        reaction_color_tuples=my_color_tuples,
        graph_attributes=my_graph_attributes,
        node_attributes=my_node_attributes,
        special_node_attributes=my_special_node_attributes,
        special_edge_attributes=my_special_edge_attributes,
    )

    t_9 = 1
    rho = 1.0e4
    mass_fractions = {
        ("n", 0, 1): 0.2,
        ("h1", 1, 1): 0.3,
        ("he4", 2, 4): 0.1,
        ("c12", 6, 12): 0.1,
        ("o16", 8, 16): 0.1,
        ("fe56", 26, 56): 0.2,
    }

    assert wn.graph.create_flow_graph(
        net,
        t_9,
        rho,
        mass_fractions,
        induced_nuc_xpath=my_induced_nuc_xpath,
        induced_reac_xpath=my_induced_reac_xpath,
        reaction_color_tuples=my_color_tuples,
        graph_attributes=my_graph_attributes,
        node_attributes=my_node_attributes,
        special_node_attributes=my_special_node_attributes,
        special_edge_attributes=my_special_edge_attributes,
    )

    my_graphs = wn.graph.create_zone_flow_graphs(
        net,
        zones,
        induced_nuc_xpath=my_induced_nuc_xpath,
        induced_reac_xpath=my_induced_reac_xpath,
        reaction_color_tuples=my_color_tuples,
        graph_attributes=my_graph_attributes,
        node_attributes=my_node_attributes,
        special_node_attributes=my_special_node_attributes,
        special_edge_attributes=my_special_edge_attributes,
    )
    for graph in my_graphs.values():
        assert graph

    my_graphs = wn.graph.create_zone_integrated_current_graphs(
        net,
        zones,
        induced_nuc_xpath=my_induced_nuc_xpath,
        induced_reac_xpath=my_induced_reac_xpath,
        reaction_color_tuples=my_color_tuples,
        graph_attributes=my_graph_attributes,
        node_attributes=my_node_attributes,
        special_node_attributes=my_special_node_attributes,
        special_edge_attributes=my_special_edge_attributes,
    )
    for graph in my_graphs.values():
        assert graph
