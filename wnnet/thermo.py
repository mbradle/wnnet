"""This module computes thermodynamic quantities for data in a nuclear
network."""

import numpy as np
import wnstatmech as ws
import wnnet.consts as wc


class Nuclei:
    """A class for computing thermodynamic quantities for a collection of
    nuclei.

    Args:
        ``nuc`` (:obj:`wnnet.nuc.Nuc`): A wnnet nuclide collection.

        ``y_cutoff`` (:obj:`float`, optional): The smallest abundance
        per nucleon to be used in thermodynamic quantity calculations.

    """

    def __init__(self, nuc, y_cutoff=1.0e-25):
        self.nuc = nuc
        self.functions = {}

        assert y_cutoff >= 0
        self.y_cut = y_cutoff

        self.functions["number density"] = self.default_number_density
        self.functions["pressure"] = self.default_pressure
        self.functions["entropy density"] = self.default_entropy_density
        self.functions["energy density"] = self.default_energy_density
        self.functions["internal energy density"] = (
            self.default_internal_energy_density
        )

    def _compute_abundances(self, mass_fractions):
        y = {}
        for key, value in mass_fractions.items():
            y_c = value / key[2]
            if y_c > self.y_cut:
                y[key[0]] = y_c
        return y

    def default_number_density(self, t9, rho, mass_fractions):
        """The default routine for computing the number density of nuclei.

        Args:
            ``t9`` (:obj:`float`):  The temperature in 10\\ :sup:`9` K at which\
            to compute the number density.

            ``rho`` (:obj:`float`):  The density in g/cc at which to compute the\
            number density.

            ``mass_fractions`` (:obj:`float`):\
                A `wnutils <https://wnutils.readthedocs.io>`_ dictionary of mass\
                fractions.

        Returns:
            A :obj:`float` giving the number density in cgs units.

        """
        assert t9 > 0
        y = self._compute_abundances(mass_fractions)

        result = 0
        for val in y.values():
            result += rho * wc.N_A * val

        return result

    def default_pressure(self, t9, rho, mass_fractions):
        """The default routine for computing the pressure of nuclei.

        Args:
            ``t9`` (:obj:`float`):  The temperature in 10\\ :sup:`9` K at which\
            to compute the pressure.

            ``rho`` (:obj:`float`):  The density in g/cc at which to compute the\
            pressure.

            ``mass_fractions`` (:obj:`float`):\
                A `wnutils <https://wnutils.readthedocs.io>`_ dictionary of mass\
                fractions.

        Returns:
            A :obj:`float` giving the pressure in cgs units.

        """
        return (
            wc.k_B
            * 1.0e9
            * t9
            * self.default_number_density(t9, rho, mass_fractions)
        )

    def default_entropy_density(self, t9, rho, mass_fractions):
        """The default routine for computing the entropy density of nuclei.

        Args:
            ``t9`` (:obj:`float`):  The temperature in 10\\ :sup:`9` K at which\
            to compute the entropy density.

            ``rho`` (:obj:`float`):  The density in g/cc at which to compute the\
            entropy density.

            ``mass_fractions`` (:obj:`float`):\
                A `wnutils <https://wnutils.readthedocs.io>`_ dictionary of mass\
                fractions.

        Returns:
            A :obj:`float` giving the entropy density in cgs units.

        """
        y = self._compute_abundances(mass_fractions)

        result = 0
        for key, value in y.items():
            result += (
                value
                * rho
                * wc.N_A
                * wc.k_B
                * (
                    5.0 / 2.0
                    - np.log(
                        value
                        / self.nuc.compute_quantum_abundance(key, t9, rho)
                    )
                )
            )

        return result

    def default_energy_density(self, t9, rho, mass_fractions):
        """The default routine for computing the energy density of nuclei,
        which includes the rest-mass energy.

        Args:
            ``t9`` (:obj:`float`):  The temperature in 10\\ :sup:`9` K at which\
            to compute the energy density (including the rest mass).

            ``rho`` (:obj:`float`):  The density in g/cc at which to compute the\
            energy density (including the rest mass).

            ``mass_fractions`` (:obj:`float`):\
                A `wnutils <https://wnutils.readthedocs.io>`_ dictionary of mass\
                fractions.

        Returns:
            A :obj:`float` giving the energy density in cgs units.

        """
        y = self._compute_abundances(mass_fractions)

        result = 0

        kt = wc.k_B * 1.0e9 * t9

        for key, value in y.items():
            result += (
                rho
                * wc.N_A
                * value
                * (
                    self.nuc.compute_atomic_mass(key) * wc.MeV_to_ergs
                    + 1.5 * kt
                )
            )

        return result

    def default_internal_energy_density(self, t9, rho, mass_fractions):
        """The default routine for computing the energy density of nuclei,
        which does not include the rest-mass energy.

        Args:
            ``t9`` (:obj:`float`):  The temperature in 10\\ :sup:`9` K at which\
            to compute the internal energy density.

            ``rho`` (:obj:`float`):  The density in g/cc at which to compute the\
            internal energy density.

            ``mass_fractions`` (:obj:`float`):\
                A `wnutils <https://wnutils.readthedocs.io>`_ dictionary of mass\
                fractions.

        Returns:
            A :obj:`float` giving the internal energy density in cgs units.

        """
        return (
            1.5
            * wc.k_B
            * 1.0e9
            * t9
            * self.default_number_density(t9, rho, mass_fractions)
        )

    def compute_quantity(self, quantity, t9, rho, mass_fractions):
        """Routine to compute the thermodynamic quantity for the nuclei.

        Args:
            ``quantity`` (:obj:`str`):  The thermodynamic quantity to compute.

            ``t9`` (:obj:`float`):  The temperature in 10\\ :sup:`9` K at which\
            to compute the quantity.

            ``rho`` (:obj:`float`):  The density in g/cc at which to compute the\
            quantity.

            ``mass_fractions`` (:obj:`float`):\
                A `wnutils <https://wnutils.readthedocs.io>`_ dictionary of mass\
                fractions.

        Returns:
            A :obj:`float` giving the thermo quantity in cgs units.

        """
        assert quantity in self.functions

        return self.functions[quantity](t9, rho, mass_fractions)

    def compute_chemical_potentials(
        self, t9, rho, mass_fractions, include_rest_mass=False
    ):
        """Routine to compute the chemical potentials of the nuclei
        in units of kT.

        Args:
            ``t9`` (:obj:`float`):  The temperature in 10\\ :sup:`9` K at which\
            to compute the chemical potentials.

            ``rho`` (:obj:`float`):  The density in g/cc at which to compute the\
            chemical potentials.

            ``mass_fractions`` (:obj:`float`):\
                A `wnutils <https://wnutils.readthedocs.io>`_ dictionary of mass\
                fractions.

            ``include_rest_mass`` (:obj:`bool`):\
                A boolean to choose whether to include the rest mass in the
                chemical potential.

        Returns:
            A :obj:`dict` with the name of each species as the key and the
            chemical potential (divide by kT) as the value.

        """
        y = self._compute_abundances(mass_fractions)

        result = {}
        for key, value in y.items():
            result[key] = np.log(
                value / self.nuc.compute_quantum_abundance(key, t9, rho)
            )
            if include_rest_mass:
                result[key] += (
                    self.nuc.compute_atomic_mass(key)
                    * wc.MeV_to_ergs
                    / (wc.k_B * 1.0e9 * t9)
                )

        return result

    def update_function(self, quantity, function):
        """Routine to update or add a thermodynamic function for the nuclei.

        Args:
            ``quantity`` (:obj:`str`):  The thermodynamic quantity.

            ``function``: A function with prototype *(t9, rho, mass_fractions)*
            to compute the thermodynamic quantity for the nuclei.

        Returns:
            On successful return, the function for the quantity has been
            updated or added.

        """
        self.functions[quantity] = function


class Thermo:
    """A class for computing thermodynamic quantities for data in a
    nuclear reaction network.

    Args:
        ``nuc`` (:obj:`wnnet.nuc.Nuc`): A wnnet nuclide collection.


    """

    def __init__(self, nuc):
        self.nuc = nuc
        self.fermions = {}
        self.bosons = {}
        self.fermions["electron"] = ws.fermion.create_electron()
        self.bosons["photon"] = ws.boson.create_photon()
        self.nuclei = Nuclei(nuc)
        self.number_density = {}

    def update_boson(self, name, boson):
        """Routine to update a boson.  If the boson does not already exist,
        it will be added.

        Args:
            ``name`` (:obj:`str`):  The name of the boson.

            ``boson``: The `wnstatmech <https://wnstatmech.readthedocs.io>`__
            boson to update.

        Returns:
            On successful return, the boson has been added or updated.

        """
        self.bosons[name] = boson

    def get_boson(self, name):
        """Routine to retrieve a boson.

        Args:
            ``name`` (:obj:`str`):  The name of the boson.

        Returns:
            The `wnstatmech <https://wnstatmech.readthedocs.io>`__ boson.

        """
        assert name in self.bosons
        return self.bosons[name]

    def remove_boson(self, name):
        """Routine to remove a boson.

        Args:
            ``name`` (:obj:`str`):  The name of the boson.

        Returns:
            On successful return, the boson has been removed.

        """
        assert name in self.bosons
        self.bosons.pop(name)

    def update_fermion(self, name, fermion):
        """Routine to update a fermion.  If the fermion does not already
        exists, it will be added.

        Args:
            ``name`` (:obj:`str`):  The name of the fermion.

            ``fermion``: The `wnstatmech <https://wnstatmech.readthedocs.io/>`__
            fermion to update.

        Returns:
            On successful return, the fermion has been added or updated.

        """
        self.fermions[name] = fermion

    def get_fermion(self, name):
        """Routine to retrieve a fermion.

        Args:
            ``name`` (:obj:`str`):  The name of the fermion.

        Returns:
            The `wnstatmech <https://wnstatmech.readthedocs.io/>`__ fermion.

        """
        assert name in self.fermions
        return self.fermions[name]

    def remove_fermion(self, name):
        """Routine to remove a fermion.

        Args:
            ``name`` (:obj:`str`):  The name of the fermion.

        Returns:
            On successful return, the fermion has been removed.

        """
        assert name in self.fermions
        self.fermions.pop(name)

    def set_number_density(self, particle, number_density):
        """Routine to update the number density for a particle.  Use this
        routine to set the number density of particles other than photons,
        electrons, and nuclei.

        Args:
            ``particle`` (:obj:`str`):  The name of the particle.

            ``number_density`` (:obj:`float`): The number density in
            cgs units.

        Returns:
            On successful return, the number density has been updated.

        """
        self.number_density[particle] = number_density

    def get_number_density(self, particle):
        """Routine to retrieve the number density for a particle.
        Use this routine for particles other than the default
        particles (photons, electrons, and nuclei).  For the default
        particles, compute the number density directly
        from the *compute_quantity* method.

        Args:
            ``particle`` (:obj:`str`):  The name of the particle.

        Returns:
            The number density of the particle in cgs units.

        """
        return self.number_density[particle]

    def compute_quantity(self, quantity, t9, rho, mass_fractions):
        """A routine to compute a thermodynamic quantity for a given set of\
        mass fractions at the input temperature and density.

        Args:
            ``quantity`` (:obj:`str`):  The thermodynamic quantity to compute.

            ``t9`` (:obj:`float`):  The temperature in 10\\ :sup:`9` K at which\
            to compute the quantity.

            ``rho`` (:obj:`float`):  The density in g/cc at which to compute the\
            quantity.

            ``mass_fractions`` (:obj:`float`):\
                A `wnutils <https://wnutils.readthedocs.io>`_ dictionary of mass\
                fractions.

        Returns:
            A :obj:`dict` with keys *electron*, *photon*, *baryon*,
            and *total* and the values being the computed quantity for
            each key.

        """

        result = {}

        temperature = 1.0e9 * t9

        # Baryons

        result["baryon"] = self.nuclei.compute_quantity(
            quantity, t9, rho, mass_fractions
        )

        # Fermions

        for name, fermion in self.fermions.items():
            if name == "electron":
                _n = 0
                for key, value in mass_fractions.items():
                    _n += rho * wc.N_A * key[1] * value / key[2]
            else:
                _n = self.get_number_density(fermion)

            alpha = fermion.compute_chemical_potential(temperature, _n)

            result[name] = fermion.compute_quantity(
                quantity, temperature, alpha
            )

        # Bosons

        for name, boson in self.bosons.items():
            if name == "photon":
                alpha = 0
            else:
                alpha = boson.compute_chemical_potential(
                    temperature, self.get_number_density(boson)
                )

            result[name] = boson.compute_quantity(quantity, temperature, alpha)

        # Total

        total = 0
        for value in result.values():
            total += value

        result["total"] = total

        return result

    def compute_quantity_in_zone(self, quantity, zone):
        """A routine to compute a thermodynamic quantity for the input zone.

        Args:
           ``quantity`` (:obj:`str`):  The thermodynamic quantity to compute.

           ``zone`` (:obj:`dict`): A `wnutils <https://wnutils.readthedocs.io>`_ *zone*.

        Returns:
            A :obj:`dict` keys *electron*, *photon*, *baryon*,
            and *total* and the values being the computed quantity for
            the zone.

        """

        return self.compute_quantity(
            quantity,
            float(zone["properties"]["t9"]),
            float(zone["properties"]["rho"]),
            zone["mass fractions"],
        )

    def compute_quantity_for_zones(self, quantity, zones):
        """A routine to compute a thermodynamic quantity for the input zones.

        Args:
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
            result[key] = self.compute_quantity_in_zone(
                quantity,
                value,
            )

        return result
