"""
chapter_2_examples.py

"""

from math import exp, log
from helpers.answer_types import SingleAnswer, MultiPartAnswer
from helpers import conversions, data
from helpers.chapter_formulae.chapter_2 import AtomicAndNuclearPhysics
from helpers.answer_writer import AnswerWriter


class Chapter2Examples:
    """Solutions of Chapter 2 example problems"""

    type_and_chapter = ("Example", 2)
    formulae = AtomicAndNuclearPhysics()

    def example_2_1(self) -> SingleAnswer:
        """A glass of water is known to contain 6.6 * 10^24 atoms of hydrogen. How
        many atoms of deuterium 2H are present?"""
        fraction = data.abundances["2H"] / 100
        n_h = 6.6e24
        n_atoms = fraction * n_h
        return SingleAnswer(ans=n_atoms, units="atoms", sig_figs=2).format()

    def example_2_2(self) -> SingleAnswer:
        """Calculate the atomic weight of naturally occuring oxygen"""
        gamma = [data.abundances["16O"], data.abundances["17O"], data.abundances["18O"]]
        m = [
            data.atomic_masses["16O"],
            data.atomic_masses["17O"],
            data.atomic_masses["18O"],
        ]
        m_avg = self.formulae.average_atomic_weight(gamma, m)
        return SingleAnswer(ans=m_avg, units=None, sig_figs=6).format()

    def example_2_3(self) -> SingleAnswer:
        """Calculate the rest-mass energy of the electron in MeV."""
        e = (
            data.constants["m_e (kg)"]
            * data.constants["c (m/s)"] ** 2
            / conversions.energy["J/MeV"]
        )
        return SingleAnswer(ans=e, units="MeV", sig_figs=4).format()

    def example_2_4(self) -> SingleAnswer:
        """Compute the energy equivalent of the atomic mass unit."""
        e = (
            self.example_2_3()[0]
            * conversions.mass["g/amu"]
            / conversions.mass["g/electron"]
        )
        return SingleAnswer(ans=e, units="MeV", sig_figs=4).format()

    def example_2_5(self) -> SingleAnswer:
        """A high-energy electron strikes a lead atom and ejects one of the K-electrons from the
        atom. What wavelength radiation is emitted when an outer electron drops into the vacancy?
        """
        ionization_energy = 88e3 * conversions.energy["J/eV"]  # J
        lambda_ = self.formulae.photon_wavelength(ionization_energy)
        return SingleAnswer(ans=lambda_, units="m", sig_figs=3).format()

    def example_2_6(self) -> MultiPartAnswer:
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
        return MultiPartAnswer(
            ans=[max_activity, t], units=["Ci", "hr"], sig_figs=3
        ).format()

    def example_2_7(self) -> str:
        """Complete the following reaction: 14N + n -> ? + 1H"""
        # The sum of atomic numbers on both sides is equal. Since Z = 1 for hydrogen, the missing
        # element is C. The total number of nucleons is supposed to be equal on both sides.
        # Therefore, Since the mass number of H is 1, the carbon isotope formed must be 14C.
        return "14N + n -> 14C + 1H"

    def example_2_8(self) -> SingleAnswer:
        """One of the reactions that occurs when 3H (tritium) is bombarded by deuterons (2H nuclei)
        is: 3H(d, n)4He. Compute the Q value of this reaction."""
        # Add 1 since the deutrons are have energy of 1 MeV
        q = (
            self.formulae.q_value_atomic_masses(
                data.atomic_masses["3H"],
                data.atomic_masses["2H"],
                data.atomic_masses["4He"],
                data.atomic_masses["n"],
            )
            + 1
        )
        return SingleAnswer(ans=q, units="MeV", sig_figs=5).format()

    def example_2_9(self) -> SingleAnswer:
        """Calculate the binding energy of the last neutron in 13C."""
        # The residual nucleus is 12C, therefore:
        m_a_1, m_a = data.atomic_masses["12C"], data.atomic_masses["13C"]
        e_s = self.formulae.separation_energy(m_a_1, m_a)
        return SingleAnswer(ans=e_s, units="MeV", sig_figs=3).format()

    def example_2_10(self) -> MultiPartAnswer:
        """Calculate the mass and binding energy of (107,47)Ag using the mass equation"""
        a, z = 107, 47
        n = a - z
        m_mev = self.formulae.mass_equation(n, a, z)
        m_amu = m_mev / conversions.energy["MeV/amu"]
        return MultiPartAnswer(
            ans=[m_mev, m_amu], units=["MeV", "amu"], sig_figs=7
        ).format()

    def example_2_11(self) -> MultiPartAnswer:
        """Calculate the most probable and average energies of air molecules in a New York City
        subway in the summer time at 38 degrees celsius"""
        t = conversions.temperature(38, units="C")  # deg. F
        e_p = self.formulae.most_probable_energy(t)
        e_p = e_p / conversions.energy["J/MeV"] * conversions.energy["eV/MeV"]
        e_avg = self.formulae.maxwellian_average_energy(t)
        e_avg = e_avg / conversions.energy["J/MeV"] * conversions.energy["eV/MeV"]
        return MultiPartAnswer(ans=[e_p, e_avg], units="eV", sig_figs=3).format()

    def example_2_12(self) -> SingleAnswer:
        """The density of sodium is 0.97 g/cm^3. Calculate its atom density"""
        rho, m = 0.97, data.atomic_masses["Na"]
        n = self.formulae.atom_density(rho, m)
        return SingleAnswer(ans=n, units="atoms/cm^3", sig_figs=3).format()

    def example_2_13(self) -> str:
        """The density of a NaCl crystal is 2.17 g/cm^3. Compute the atom densities of Na and Cl."""
        # Since there is one Na atom and one Cl atom in a pseudo molecule, add the weights:
        m = data.atomic_masses["Na"] + data.atomic_masses["Cl"]
        rho = 2.17
        n = self.formulae.atom_density(rho, m)
        return f"{n:.3g} atoms/cm^3 for both since there is one atom of each per pseudomolecule"

    def example_2_14(self) -> MultiPartAnswer:
        """For water of normal (unit) density compute:
        (a) the number of H2O molecules per cm^3,
        (b) the atom densities of hydrogen and oxygen
        (c) the atom density of 2H"""
        # a)
        m = 2 * data.atomic_masses["H"] + data.atomic_masses["O"]
        rho = 1
        n_h2o = self.formulae.atom_density(rho, m)
        # b) Multiply result of part a) by 2 for hydrogen, oxygen stays the same.
        n_h = n_h2o * 2
        n_o = n_h2o
        # c)
        gamma = data.abundances["2H"]
        n_2h = gamma * n_h / 100
        return MultiPartAnswer(
            ans=[n_h2o, n_h, n_o, n_2h],
            units=["molecules/cm^3", "atoms/cm^3", "atoms/cm^3", "atoms/cm^3"],
            sig_figs=4,
        ).format()

    def example_2_15(self) -> MultiPartAnswer:
        """A certain nuclear reactor is fueled with 1,500 kg of uranium rods enriched to 20 w/0
        235U. The remainder is 238U. The density of the uranium is 19.1 g/cm^3
        (a) How much 235U is in the reactor?
        (b) What are the atom densities of 235U and 238U in the rods?"""
        m_total, rho = 1500, 19.1
        m_235u, m_238u = data.atomic_masses["235U"], data.atomic_masses["238U"]
        # a) Enrichment to 20% w/o means 20% of the uranium in the reactor is 235U
        wp_235u, wp_238u = 20, 80
        mass_235u = 0.2 * m_total
        # b)
        n_235u = self.formulae.atom_density_component(wp_235u, rho, m_235u)
        n_238u = self.formulae.atom_density_component(wp_238u, rho, m_238u)
        return MultiPartAnswer(
            ans=[mass_235u, n_235u, n_238u],
            units=["kg", "atoms/cm^3", "atoms/cm^3"],
            sig_figs=3,
        ).format()

    def example_2_16(self) -> SingleAnswer:
        """The fuel for a reactor consists of pellets of uranium dioxide (UO2), which have a density
        of 10.5 g/cm^3. If the uranium is enriched to 30 w/o in 235U, what is the atom density of
        the 235U in the fuel?"""
        rho = 10.5
        wp_235u, wp_238u = 30, 70
        m_235u, m_238u = data.atomic_masses["235U"], data.atomic_masses["238U"]
        m_u = self.formulae.atom_density_weight_percent(
            [wp_235u, wp_238u], [m_235u, m_238u]
        )
        m_o = data.atomic_masses["O"]
        wp_u = self.formulae.weight_percent(m_u, 1, m_o, 2)
        rho_avg_u = wp_u * rho / 100
        rho_235u = wp_235u * rho_avg_u / 100
        n_235u = self.formulae.atom_density(rho_235u, m_235u)
        return SingleAnswer(ans=n_235u, units="atoms/cm^3", sig_figs=3).format()


if __name__ == "__main__":
    ch2e = Chapter2Examples()
    writer = AnswerWriter(obj=ch2e, num_questions=16, path="./answers/examples.txt")
    writer.write_answers()
