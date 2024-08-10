import requests, io
import wnnet as wn


def get_net():
    return wn.net.Net(
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
