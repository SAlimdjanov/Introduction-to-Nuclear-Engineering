"""
chapter_2.py

TODO:
Method q_value_atomic_masses: Change input args to lists of reactants and products if necessary

"""

from math import sqrt, exp, log, pi
from typing import List
from .. import conversions, data


class AtomicAndNuclearPhysics:
    """Chapter 2"""

    def average_atomic_weight(self, gamma: List[float], m: List[float]) -> float:
        """Computes the average atomic weight of a mixture of isotopes, denoted by M. Requires lists
        of abundances and atomic masses of each isotope in the correct order. len(gamma) must be
        equal to len(m)."""
        n = len(gamma)
        total = sum(gamma[i] * m[i] for i in range(n))
        return total / 100

    def abundance_natural(
        self, gamma_i: List[float], m_i: List[float], gamma: float, m: float
    ) -> float:
        """Computes the abundance of a nuclide in the naturally occuring form of its element, given
        the masses and abundances of all possibilities and the abundance and mass of the target
        nuclide"""
        return gamma * m / sum(gamma_i[j] * m_i[j] for j in range(len(gamma_i)))

    def nuclear_radius(self, a: float) -> float:
        """Computes a nucleus's radius in fm given its atomic mass number"""
        return 1.25 * pow(a, 1 / 3)

    def rest_energy(self, m_0: float) -> float:
        """Computes the rest energy of a particle given its rest mass"""
        return m_0 * data.constants["c (m/s)"] ** 2

    def relative_mass(self, m_0: float, v: float) -> float:
        """Calculates the relative mass of a moving particle"""
        return m_0 / sqrt(1 - v**2 / data.constants["c (m/s)"] ** 2)

    def total_energy(self, m: float) -> float:
        """Calculates the total energy of a particle given its relative mass"""
        return m * data.constants["c (m/s)"] ** 2

    def kinetic_energy_rel(self, m_0: float, v: float) -> float:
        """Computes a particle's kinetic energy given its rest mass and speed"""
        return (self.relative_mass(m_0, v) - m_0) * data.constants["c (m/s)"] ** 2

    def kinetic_energy(self, m_0: float, v: float) -> float:
        """Calculates a particle's non-relativistic kinetic energy"""
        return 0.5 * m_0 * v**2

    def particle_speed(self, e_t: float, e_r: float) -> float:
        """Computes the speed of any particle, relativistic or nonrelativistic, given its total
        energy and rest energy"""
        c = data.constants["c (m/s)"]
        return c * sqrt(1 - e_r**2 / e_t**2)

    def neutron_speed(self, e: float) -> float:
        """Computes a neutron's speed in cm/s given its kinetic energy in eV"""
        return 1.383e6 * sqrt(e)

    def photon_energy(self, nu: float) -> float:
        """Obtains a photon's energy given its frequency. Second arg selects the units of Planck's
        constant. Can choose 'J-s' or 'eV-s'"""
        h = data.constants["h (J-s)"]
        return h * nu

    def photon_wavelength(self, e: float) -> float:
        """Calculates a photon's wavelength given its energy"""
        h, c = data.constants["h (J-s)"], data.constants["c (m/s)"]
        return h * c / e

    def wavelength(self, p: float) -> float:
        """Computes a particle's wavelength given its momentum"""
        h = data.constants["h (J-s)"]
        return h / p

    def wavelength_nonrel(self, m_0: float, e: float) -> float:
        """Calculates a particle's non-relativistic wavelength given its rest mass and kinetic
        energy"""
        h = data.constants["h (J-s)"]
        return h / sqrt(2 * m_0 * e)

    def wavelength_w_compton(self, e_t: float, e_r: float) -> float:
        """Computes a particle's relativistic wavelength given its total energy and rest energy,
        in relation to the Compton wavelength"""
        m_e, h, c = (
            data.constants["m_e (kg)"],
            data.constants["h (J-s)"],
            data.constants["c (m/s)"],
        )
        lambda_c = h / (m_e * c)
        return lambda_c * m_e * c**2 / sqrt(e_t**2 - e_r**2)

    def neutron_wavelength(self, e: float) -> float:
        """Returns a neutron's wavelength in cm given its kinetic energy in eV"""
        return 2.860e-9 / sqrt(e)

    def momentum_rel(self, e_t: float, e_r: float) -> float:
        """Computes the relativistic momentum of a particle given its total energy and rest
        energy"""
        return sqrt(e_t**2 - e_r**2) / data.constants["c (m/s)"]

    def wavelength_rel(self, e_t: float, e_r: float) -> float:
        """Computes the relativistic wavelength of a particle given its total energy and rest
        energy"""
        h = data.constants["h (J-s)"]
        return h * data.constants["c (m/s)"] / sqrt(e_t - e_r)

    def decay_constant(self, half_life: float) -> float:
        """Obtains decay constant given an element's half-life"""
        return log(2) / half_life

    def undecayed_atoms(self, n_0: float, lambda_: float, t: float) -> float:
        """Calculates the number of undecayed atoms present in a sample"""
        return n_0 * exp(-lambda_ * t)

    def activity(self, alpha_0: float, lambda_: float, t: float) -> float:
        """Computes the activity of the sample given its initial activity, decay constant, and the
        point in time during decay"""
        return alpha_0 * exp(-lambda_ * t)

    def activity_alt(self, alpha_0: float, half_life: float, t: float) -> float:
        """Calculates the activity at t if lambda is not known"""
        return alpha_0 * exp(-log(2) * t / half_life)

    def activity_atom_density(self, lambda_: float, m: float, gamma: float) -> float:
        """Computes the activity of an element given its decay constant, atomic weights, and
        abundance"""
        n_a = data.constants["N_A (1/(g-mol))"]
        return lambda_ * n_a * gamma / m

    def half_life(self, lambda_: float) -> float:
        """Returns the half-life of an element given its decay constant"""
        return log(2) / lambda_

    def mean_life(self, lambda_=None, half_life=None):
        """Computes the mean-life of a sample given its decay constant or half-life"""
        if lambda_ is None and half_life is None:
            raise ValueError("Half-life or decay constant is required")
        return 1 / lambda_ if lambda_ else half_life / log(2)

    def undecayed_atoms_reactor(
        self, n_0: float, r: float, lambda_: float, t: float
    ) -> float:
        """For problems where atoms are produced in a reactor at rate R (atoms/s)"""
        return n_0 * exp(-lambda_ * t) + r * (1 - exp(-lambda_ * t)) / lambda_

    def activity_reactor_atoms(
        self, n_0: float, r: float, lambda_: float, t: float
    ) -> float:
        """Reactor activity with the number of atoms initial present"""
        return lambda_ * self.undecayed_atoms_reactor(n_0, r, lambda_, t)

    def activity_reactor(
        self, alpha_0: float, r: float, lambda_: float, t: float
    ) -> float:
        """Reactor activity with the initial activity value, alpha_0"""
        return alpha_0 * exp(-lambda_ * t) + r * (1 - exp(-lambda_ * t))

    def number_of_atoms_decay_chain(
        self, n_a0: float, n_b0: float, lambda_a: float, lambda_b: float, t: float
    ) -> float:
        """Calculates the number of atoms of element B after decaying from element A at time t"""
        return n_b0 * exp(-lambda_b * t) + n_a0 * lambda_a / (lambda_b - lambda_a) * (
            exp(-lambda_a * t) - exp(-lambda_b * t)
        )

    def activity_decay_chain(
        self,
        alpha_a0: float,
        alpha_b0: float,
        lambda_a: float,
        lambda_b: float,
        t: float,
    ) -> float:
        """Calculates the activity of element B after decaying from element A at time t"""
        return alpha_b0 * exp(-lambda_b * t) + alpha_a0 * lambda_b / (
            lambda_b - lambda_a
        ) * (exp(-lambda_a * t) - exp(-lambda_b * t))

    def q_value(
        self,
        m_a: float,
        z_a: int,
        m_b: float,
        z_b: int,
        m_c: float,
        z_c: int,
        m_d: float,
        z_d: int,
    ) -> float:
        """Computes the q value with atomic masses and atomic numbers in MeV"""
        m_e = data.constants["m_e (MeV)"]
        return (
            ((m_a + z_a * m_e) + (m_b + z_b * m_e))
            - ((m_c + z_c * m_e) + (m_d + z_d * m_e))
        ) * conversions.energy["MeV/amu"]

    def q_value_atomic_masses(
        self, m_a: float, m_b: float, m_c: float, m_d: float
    ) -> float:
        """Computes the q value with atomic masses in MeV"""
        return ((m_a + m_b) - (m_c + m_d)) * conversions.energy["MeV/amu"]

    def mass_defect(self, z: int, n: int, m: int) -> float:
        """Computes the mass defect of a nucleus"""
        m_h, m_n = data.atomic_masses["1H"], data.atomic_masses["n"]
        return z * m_h + n * m_n - m

    def q_value_binding(
        self, be_a: float, be_b: float, be_c: float, be_d: float
    ) -> float:
        """Computes the Q value with binding energies"""
        return (be_c + be_d) - (be_a - be_b)

    def q_value_mass_defects(
        self, d_a: float, d_b: float, d_c: float, d_d: float
    ) -> float:
        """Obtains the Q value of a reaction given the mass defects of the reactants and products
        in MeV"""
        return (d_a + d_b) - (d_c + d_d)

    def binding_energy(self, m: float, z: int, n: int) -> float:
        """Computes the binding energy of a nuclide given its mass in amu, atomic number, and
        neutron number"""
        m_n, m_1h = data.atomic_masses["n"], data.atomic_masses["1H"]
        return (z * m_1h + n * m_n - m) * conversions.energy["MeV/amu"]

    def separation_energy(self, m_a_1: float, m_a: float) -> float:
        """Calculates the separation energy of the last neutron in nucleus ^{A}Z with the mass of
        the residual nucleus ^{A-1}Z in MeV"""
        m_n = data.atomic_masses["n"]
        return (m_n + m_a_1 - m_a) * conversions.energy["MeV/amu"]

    def mass_equation(self, n: int, a: int, z: int) -> float:
        """Computes the mass of a nuclide with the liquid drop model mass equation in MeV"""
        m_n, m_p, alpha, beta, gamma, zeta, delta = (
            data.constants["m_n (MeV)"],
            data.constants["m_p (MeV)"],
            data.constants["alpha (MeV)"],
            data.constants["beta (MeV)"],
            data.constants["gamma (MeV)"],
            data.constants["zeta (MeV)"],
            data.constants["delta (MeV)"],
        )
        if (n % 2 != 0 and z % 2 == 0) or (n % 2 == 0 and z % 2 != 0):
            delta_term = 0
        elif n == 1 and z == 1:
            delta_term = 0
        elif n % 2 == 0 and z % 2 == 0:
            delta_term = -delta
        else:
            delta_term = delta
        return (
            n * m_n
            + z * m_p
            - alpha * a
            + beta * a ** (2 / 3)
            + gamma * z**2 / a ** (1 / 3)
            + zeta * ((a - 2 * z) ** 2) / a
            + delta_term
        )

    def binding_energy_mass_eqn(self, n: int, a: int, z: int) -> float:
        """Computes an approximation of the total binding energy using the last five terms of the
        mass equation, given the neutron number, atomic mass number, and atomic number
        """
        alpha, beta, gamma, zeta, delta = (
            data.constants["alpha (MeV)"],
            data.constants["beta (MeV)"],
            data.constants["gamma (MeV)"],
            data.constants["zeta (MeV)"],
            data.constants["delta (MeV)"],
        )
        if (n % 2 != 0 and z % 2 == 0) or (n % 2 == 0 and z % 2 != 0):
            delta_term = 0
        elif n == 1 and z == 1:
            delta_term = 0
        elif n % 2 == 0 and z % 2 == 0:
            delta_term = -delta
        else:
            delta_term = delta
        return -(
            -alpha * a
            + beta * a ** (2 / 3)
            + gamma * z**2 / a ** (1 / 3)
            + zeta * ((a - 2 * z) ** 2) / a
            + delta_term
        )

    def maxwellian_energy_dist(self, n: float, t: float, e: float) -> float:
        """Computes the density of particles per unit energy, N(E) given the number of particles
        N, temperature T, and energy E."""
        k = data.constants["k_B (J/K)"]
        return 2 * pi * n / pow(pi * k * t, 3 / 2) * sqrt(e) * exp(-e / (k * t))

    def most_probable_energy(self, t: float) -> float:
        """Returns the most probable energy in a Maxwellian distribution at a given temperature in
        degrees Kelvin"""
        k = data.constants["k_B (J/K)"]
        return k * t / 2

    def maxwellian_average_energy(self, t: float) -> float:
        """Calculates the average energy across a Maxwellian distribution function at a given
        temperature"""
        return 3 * self.most_probable_energy(t)

    def ideal_gas_pressure(self, n: float, t: float) -> float:
        """Computes the pressure of a gas using the ideal gas law, given its density (atoms/m^3)
        and temperature"""
        k = data.constants["k_B (J/K)"]
        return n * k * t

    def atom_density(self, rho: float, m: float) -> float:
        """Returns the atom density of an element given its physical density and gram atomic
        weight"""
        n_a = data.constants["N_A (1/(g-mol))"]
        return rho * n_a / m

    def atom_density_isotope(self, gamma: float, rho: float, m: float) -> float:
        """Calculates the atom density of an isotope given its abundance, physical density and
        atomic weight"""
        n_a = data.constants["N_A (1/(g-mol))"]
        return gamma * rho * n_a / (m * 100)

    def average_density_component(self, w: float, rho: float) -> float:
        """Computes the average density of a component given its weight percent and the physical
        density of the mixture it is contained in"""
        return w * rho / 100

    def atom_density_component(self, w: float, rho: float, m: float) -> float:
        """Calculates the atom density of a component given its weight percent, physical
        density and atomic weight"""
        n_a = data.constants["N_A (1/(g-mol))"]
        return w * rho * n_a / (m * 100)

    def weight_percent(self, x: float, m: int, y: float, n: int) -> float:
        """Calculates the weight percent of element X in a compound X_m Y_n, where args x and y
        are the atomic weights of X and Y. m and n are the number of x and y atoms in the
        compound."""
        return 100 * x * m / (x * m + y * n)

    def atom_density_w_densities(self, n: List[float]) -> float:
        """Computes the atom density of an element given the densities of its isotopes."""
        return sum(n)

    def atom_density_weight_percent(self, w: List[float], m: List[float]) -> float:
        """Returns the atom density of an element given the weight percents and atomic weights of
        its isotopes"""
        n = len(w)
        sum_ = sum(w[i] / m[i] for i in range(n)) / 100
        return 1 / sum_
