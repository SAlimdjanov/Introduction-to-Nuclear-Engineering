"""
chapter_2_examples.py

"""

from math import exp, log
from typing import List, Tuple
from helpers.data import constants, abundances, atomic_weights, neutral_atomic_masses
from helpers.conversions import conversion_factors, convert_temperature
from helpers.chapter_formulae.chapter_2 import AtomicAndNuclearPhysics
from helpers.answer_writer import AnswerWriter


class Chapter2Examples:
    """Solutions of Chapter 2 example problems"""

    formulae = AtomicAndNuclearPhysics()

    def __init__(self) -> None:
        self.type_and_chapter = ("Example", 2)

    def example_2_1(self) -> Tuple[float, str, int]:
        """A glass of water is known to contain 6.6 * 10^24 atoms of hydrogen. How
        many atoms of deuterium 2H are present?"""
        fraction = abundances["2H"] / 100
        number_of_atoms = 6.6e24
        return fraction * number_of_atoms, "atoms", 2

    def example_2_2(self) -> Tuple[float, str, int]:
        """Calculate the atomic weight of naturally occuring oxygen"""
        gamma = [abundances["16O"], abundances["17O"], abundances["18O"]]
        m = [
            atomic_weights["16O"],
            atomic_weights["17O"],
            atomic_weights["18O"],
        ]
        return self.formulae.average_atomic_weight(gamma, m), "", 6

    def example_2_3(self) -> Tuple[float, str, int]:
        """Calculate the rest-mass energy of the electron in MeV."""
        energy = constants["m_e (kg)"] * constants["c (m/s)"] ** 2
        return energy / conversion_factors["J/MeV"], "MeV", 4

    def example_2_4(self) -> Tuple[float, str, int]:
        """Compute the energy equivalent of the atomic mass unit."""
        energy = self.example_2_3()[0]
        return (
            energy * conversion_factors["g/amu"] / conversion_factors["g/electron"],
            "MeV",
            4,
        )

    def example_2_5(self) -> Tuple[float, str, int]:
        """A high-energy electron strikes a lead atom and ejects one of the K-electrons from the
        atom. What wavelength radiation is emitted when an outer electron drops into the vacancy?
        """
        ionization_energy = 88e3  # 88 keV
        return 1.240e-6 / ionization_energy, "m", 3

    def example_2_6(self) -> List[Tuple[float, str, int]]:
        """Gold-198 (half-life 64.8 hr) can be produced by bombarding stable 197Au with neutrons in
        a nuclear reactor. Suppose that a 197Au foil weighing 0.1 g is placed in a certain reactor
        for 12 hrs and that its activity is 0.90 Ci when removed.
        (a) What is the theoretical maximum activity due to 198Au in the foil?
        (b) How long does it take for the activity to reach 80 percent of the maximum?
        """
        lambda_ = self.formulae.decay_constant(half_life=64.8)  # Units: atoms/hr
        # a) Rearrange activity equation for R given alpha_0 = 0
        max_activity = 0.9 / (1 - exp(-lambda_ * 12))
        # b) Now that R is known, rearrange activity equation for t
        t = -log(0.2) / lambda_
        return [(max_activity, "Ci", 3), (t, "hr", "2")]

    def example_2_7(self) -> str:
        """Complete the following reaction: 14N + n -> ? + 1H"""
        # The sum of atomic numbers on both sides is equal. Since Z = 1 for hydrogen, the missing
        # element is C. The total number of nucleons is supposed to be equal on both sides.
        # Therefore, Since the mass number of H is 1, the carbon isotope formed must be 14C.
        return "14N + n -> 14C + 1H"

    def example_2_8(self) -> Tuple[float, str, int]:
        """One of the reactions that occurs when 3H (tritium) is bombarded by deuterons (2H nuclei)
        is: 3H(d, n)4He. Compute the Q value of this reaction."""
        # Add 1 since the deutrons are have energy of 1 MeV
        q = (
            self.formulae.q_value_atomic_masses(
                neutral_atomic_masses["3H"],
                neutral_atomic_masses["2H"],
                neutral_atomic_masses["4He"],
                neutral_atomic_masses["n"],
            )
            + 1
        )
        return q, "MeV", 5

    def example_2_9(self) -> Tuple[float, str, int]:
        """Calculate the binding energy of the last neutron in 13C."""
        # The residual nucleus is 12C, therefore:
        m_a_1, m_a = neutral_atomic_masses["12C"], neutral_atomic_masses["13C"]
        e_s = self.formulae.separation_energy(m_a_1, m_a)
        return e_s, "MeV", 3

    def example_2_10(self) -> List[Tuple[float, str, int]]:
        """Calculate the mass and binding energy of (107,47)Ag using the mass equation"""
        a, z = 107, 47
        n = a - z
        mass_mev = self.formulae.mass_equation(n, a, z)
        mass_amu = mass_mev / conversion_factors["MeV/amu"]
        return [(mass_mev, "MeV", 7), (mass_amu, "amu", 7)]

    def example_2_11(self) -> List[Tuple[float, str, int]]:
        """Calculate the most probable and average energies of air molecules in a New York City
        subway in the summer time at 38 degrees celsius"""
        t = convert_temperature(38, units="C")
        e_p = self.formulae.most_probable_energy(t)
        e_p = e_p / conversion_factors["J/MeV"] * conversion_factors["eV/MeV"]
        e_avg = self.formulae.maxwellian_average_energy(t)
        e_avg = e_avg / conversion_factors["J/MeV"] * conversion_factors["eV/MeV"]
        return [(e_p, "eV", 3), (e_avg, "eV", 3)]

    def example_2_12(self) -> Tuple[float, str, int]:
        """The density of sodium is 0.97 g/cm^3. Calculate its atom density"""
        rho, m = 0.97, atomic_weights["Na"]
        n = self.formulae.atom_density(rho, m)
        return n, "atoms/cm^-3", 3

    def example_2_13(self) -> str:
        """The density of a NaCl crystal is 2.17 g/cm^3. Compute the atom densities of Na and Cl."""
        # Since there is one Na atom and one Cl atom in a pseudo molecule, add the weights:
        m = atomic_weights["Na"] + atomic_weights["Cl"]
        rho = 2.17
        n = self.formulae.atom_density(rho, m)
        return f"{n:.3g} atoms/cm^3 for both since there is one atom of each per pseudomolecule"

    def example_2_14(self) -> List[Tuple[float, str, int]]:
        """For water of normal (unit) density compute:
        (a) the number of H2O molecules per cm^3,
        (b) the atom densities of hydrogen and oxygen
        (c) the atom density of 2H"""
        # a)
        m = 2 * atomic_weights["H"] + atomic_weights["O"]
        rho = 1
        n_h2o = self.formulae.atom_density(rho, m)
        # b) Multiply result of part a) by 2 for hydrogen, oxygen stays the same.
        n_h = n_h2o * 2
        n_o = n_h2o
        # c)
        gamma = abundances["2H"]
        n_2h = gamma * n_h / 100
        return [
            (n_h2o, "molecules/cm^3", 4),
            (f"N(H) = {n_h:.4g} and N(O) = {n_o:.4g} atoms/cm^3"),
            (n_2h, "atoms/cm^3", 5),
        ]

    def example_2_15(self) -> str:
        """A certain nuclear reactor is fueled with 1,500 kg of uranium rods enriched to 20 w/0
        235U. The remainder is 238U. The density of the uranium is 19.1 g/cm^3
        (a) How much 235U is in the reactor?
        (b) What are the atom densities of 235U and 238U in the rods?"""
        m_total, rho = 1500, 19.1
        m_235u, m_238u = atomic_weights["235U"], atomic_weights["238U"]
        # a) Enrichment to 20% w/o means 20% of the uranium in the reactor is 235U
        wp_235u, wp_238u = 20, 80
        mass_235u = 0.2 * m_total
        # b)
        n_235u = self.formulae.atom_density_component(wp_235u, rho, m_235u)
        n_238u = self.formulae.atom_density_component(wp_238u, rho, m_238u)
        return [
            (mass_235u, "kg", 4),
            (f"N(235U) = {n_235u:.3g} and N(238U) = {n_238u:.3g} atoms/cm^-3"),
        ]

    def example_2_16(self) -> str:
        """The fuel for a reactor consists of pellets of uranium dioxide (UO2), which have a density
        of 10.5 g/cm^3. If the uranium is enriched to 30 w/o in 235U, what is the atom density of
        the 235U in the fuel?"""
        rho = 10.5
        wp_235u, wp_238u = 30, 70
        m_235u, m_238u = atomic_weights["235U"], atomic_weights["238U"]
        m_u = self.formulae.atom_density_weight_percent(
            [wp_235u, wp_238u], [m_235u, m_238u]
        )
        m_o = atomic_weights["O"]
        wp_u = self.formulae.weight_percent(m_u, 1, m_o, 2)
        rho_avg_u = wp_u * rho / 100
        rho_235u = wp_235u * rho_avg_u / 100
        n_235u = self.formulae.atom_density(rho_235u, m_235u)
        return n_235u, "atoms/cm^3", 3


if __name__ == "__main__":
    ch2e = Chapter2Examples()
    ans = AnswerWriter(ch2e, 16, "./answers/examples.txt")
    ans.write_answers()
