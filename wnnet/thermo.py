"""This module computes thermodynamic quantities for a network."""

import numpy as np
import wnstatmech as ws
import wnnet.consts as wc


def _compute_baryon_quantity(nuc, quantity, t_9, rho, mass_fractions):

    result = 0

    kt = wc.k_B * 1.0e9 * t_9

    for key, value in mass_fractions.items():
        y = value / key[2]

        if y > 1.e100:

            n_den = rho * wc.N_A * y

            if quantity == "pressure":
                result += n_den * kt
            elif quantity == "energy density":
                result += n_den * (
                    nuc.compute_atomic_mass(key[0]) * wc.MeV_to_ergs
                    + 1.5 * kt
                )
            elif quantity == "internal energy density":
                result += n_den * 1.5 * kt
            elif quantity == "entropy density":
                result += (
                    n_den
                    * wc.k_B
                    * (
                        5.0 / 2.0
                        - np.log(
                            y / nuc.compute_quantum_abundance(key[0], t_9, rho)
                        )
                    )
                )

    return result


def compute_thermo_quantity(nuc, quantity, t_9, rho, mass_fractions):
    """A routine to compute a thermodynamic quantity for a given set of\
    mass fractions at the input temperature and density.

    Args:
        ``nuc``: A wnnet nuclide collection.

        ``quantity`` (:obj:`str`):  The thermodynamic quantity to compute.

        ``t_9`` (:obj:`float`):  The temperature in 10\\ :sup:`9` K at which\
          to compute the flows.

        ``rho`` (:obj:`float`):  The density in g/cc at which to compute the\
          flows.

        ``mass_fractions`` (:obj:`float`):\
        A `wnutils <https://wnutils.readthedocs.io>`_ dictionary of mass\
        fractions.

    Returns:
        A :obj:`dict` with keys *electron*, *photon*, *baryon*,
        and *total* and the values being the computed quantity for
        each key.

    """

    assert quantity in (
        "pressure",
        "energy density",
        "internal energy density",
        "entropy density",
    )

    result = {}

    temperature = 1.0e9 * t_9

    # Electrons

    y_e = 0
    for key, value in mass_fractions.items():
        y_e += key[1] * value / key[2]

    n_e = rho * wc.N_A * y_e

    electron = ws.fermion.create_electron()

    alpha = electron.compute_chemical_potential(temperature, n_e)

    result["electron"] = electron.compute_quantity(
        quantity, temperature, alpha
    )

    # Photons

    photon = ws.boson.create_photon()

    result["photon"] = photon.compute_quantity(quantity, temperature, 0)

    # Baryons

    result["baryon"] = _compute_baryon_quantity(
        nuc, quantity, t_9, rho, mass_fractions
    )

    # Total

    total = 0
    for value in result.values():
        total += value

    result["total"] = total

    return result


def compute_thermo_quantity_in_zone(nuc, quantity, zone):
    """A routine to compute a thermodynamic quantity for the input zone.

    Args:
        ``nuc``: A wnnet nuclide collection.

        ``quantity`` (:obj:`str`):  The thermodynamic quantity to compute.

        ``zone`` (:obj:`dict`): A `wnutils <https://wnutils.readthedocs.io>`_ *zone*.

    Returns:
        A :obj:`dict` keys *electron*, *photon*, *baryon*,
        and *total* and the values being the computed quantity for
        the zone.

    """

    return compute_thermo_quantity(
        nuc,
        quantity,
        float(zone["properties"]["t9"]),
        float(zone["properties"]["rho"]),
        zone["mass fractions"],
    )


def compute_thermo_quantity_for_zones(nuc, quantity, zones):
    """A routine to compute a thermodynamic quantity for the input zone.

    Args:
        ``nuc``: A wnnet nuclide collection.

        ``quantity`` (:obj:`str`):  The thermodynamic quantity to compute.

        ``zones`` (:obj:`dict`): A dictionary of
         `wnutils <https://wnutils.readthedocs.io>`_ *zone data*.

    Returns:
        A :obj:`dict` with the key being the zone label and the value
        being another dictionary with keys *electron*, *photon*, *baryon*,
        and *total* and the values being the computed quantity for
        each zone.

    """

    result = {}

    for key, value in zones.items():
        result[key] = compute_thermo_quantity_in_zone(
            nuc,
            quantity,
            value,
        )

    return result
