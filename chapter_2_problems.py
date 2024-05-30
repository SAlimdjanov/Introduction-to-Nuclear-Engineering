"""
chapter_2_problems.py

"""

from math import log, sqrt, pi
from typing import List, Tuple
from helpers.data import constants, abundances, atomic_masses, half_lives
from helpers.conversions import conversion_factors
from helpers.chapter_formulae.chapter_2 import AtomicAndNuclearPhysics
from helpers.answer_writer import AnswerWriter


class Chapter2Problems:
    """Solutions to Chapter 2 Problems"""

    type_and_chapter = ("Problem", 2)
    formulae = AtomicAndNuclearPhysics()

    def problem_2_1(self) -> List[str]:
        """How many neutrons and protons are there in the nuclei of the following atoms:
        (a) 7Li, (b) 24Mg, (c) 135Xe, (d) 209Bi, (e) 222Rn"""
        return [
            "n = 4 and p = 3",
            "n = 12 and p = 12",
            "n = 81 and p = 54",
            "n = 126 and p = 83",
            "n = 136 and p = 86",
        ]

    def problem_2_2(self) -> str:
        """The atomic weight of 59Co is 58.93319. How many times heavier is 59Co than 12C?"""
        ratio = atomic_masses["59Co"] / atomic_masses["12C"]
        return f"59Co is {ratio:.3g} times heavier than 12C"

    def problem_2_3(self) -> Tuple[float, str, int]:
        """How many atoms are there in 10 g of 12C?"""
        mass_12c = 10
        mol_12c = mass_12c / atomic_masses["12C"] * abundances["12C"] / 100
        n_a = constants["N_A (1/(g-mol))"]
        atoms = mol_12c * n_a
        return atoms, "atoms", 4

    def problem_2_4(self) -> List[Tuple[float, str, int]]:
        """Compute the molecular weights of:
        (a) H_2 gas, (b) H_2 O, (c) H_2 O_2"""
        gamma_1h, gamma_2h = abundances["1H"], abundances["2H"]
        m_1h, m_2h = atomic_masses["1H"], atomic_masses["2H"]
        # a) hydrogen gas
        m_h2 = (
            self.formulae.average_atomic_weight([gamma_1h, gamma_2h], [m_1h, m_2h]) * 2
        )
        # b) Water
        gamma_i = [abundances["16O"], abundances["17O"], abundances["18O"]]
        m_i = [atomic_masses["16O"], atomic_masses["17O"], atomic_masses["18O"]]
        m_o = self.formulae.average_atomic_weight(gamma_i, m_i)
        m_h2o = m_h2 + m_o
        # c) Hydrogen peroxide
        m_h2o2 = m_h2 + m_o * 2
        return [(m_h2, "amu", 8), (m_h2o, "amu", 8), (m_h2o2, "amu", 8)]

    def problem_2_5(self) -> List[Tuple[float, str, int]]:
        """When H_2 gas is formed from naturally occurring hydrogen, what percentages of the
        molecules have molecular weights of approximately 2, 3, and 4?"""
        # Four possible combinations involving hydrogen and deuterium: 1H 1H, 2H 1H, 1H 2H, 2H 2H
        gamma_1h, gamma_2h = abundances["1H"], abundances["2H"]
        mass_2 = ((gamma_1h / 100) ** 2) * 100
        mass_3 = gamma_1h * gamma_2h / 100 * 2
        mass_4 = ((gamma_2h / 100) ** 2) * 100
        return [(mass_2, "%", 4), (mass_3, "%", 4), (mass_4, "%", 4)]

    def problem_2_6(self) -> Tuple[float, str, int]:
        """Natural uranium is composed of three isotopes: 234U, 235U, and 238U. Their abundances
        and atomic weights are given in the following table. Compute the atomic weight of natural
        uranium."""
        gamma_i = [abundances["234U"], abundances["235U"], abundances["238U"]]
        m_i = [atomic_masses["234U"], atomic_masses["235U"], atomic_masses["238U"]]
        m = self.formulae.average_atomic_weight(gamma_i, m_i)
        return m, "amu", 8

    def problem_2_7(self) -> List[Tuple[float, str, int]]:
        """A beaker contains 50 g of naturally occuring water
        (a) How many moles of water are present?
        (b) How many atoms are hydrogen atoms?
        (c) How many atoms are deuterium atoms?"""
        mass_h2o = 50
        # a) Obtain weight from problem 2.4 b)
        moles = mass_h2o / self.problem_2_4()[1][0]
        # b) There are two hydrogen atoms per water molecule:
        hydrogen_atoms = 2 * moles * constants["N_A (1/(g-mol))"]
        # c) Multiply by the deuterium abundance fraction
        deuterium_atoms = abundances["2H"] * hydrogen_atoms / 100
        return [
            (moles, "mol", 5),
            (hydrogen_atoms, "atoms", 5),
            (deuterium_atoms, "atoms", 5),
        ]

    def problem_2_8(self) -> Tuple[float, str, int]:
        """The glass in Example 2.1 has an inside diameter of 7.5 cm. How high does the water stand
        in the glass?"""
        h_atoms = 6.6e24
        # If there are x hydrogen atoms in a glass of water, there are x / 2 oxygen atoms
        o_atoms = h_atoms / 2
        molecules = (h_atoms + o_atoms) / 3
        mol_h2o = molecules / constants["N_A (1/(g-mol))"]
        m_h2o = self.problem_2_4()[1][0]  # Molecular mass of water
        mass_h2o = m_h2o * mol_h2o
        # The density is given by m * V. Substituting V with the volume of a cylinder and solve for
        # the height:
        rho_h2o, r = 1, 7.5 / 2  # g/cm^3, cm
        h = mass_h2o / (rho_h2o * pi * r**2)
        return h, "cm", 4

    def problem_2_9(self) -> Tuple[float, str, int]:
        """Calculate the mass of a proton in amu"""
        m_p = (
            constants["m_p (kg)"]
            * conversion_factors["g/kg"]
            / conversion_factors["g/amu"]
        )
        return m_p, "amu", 5

    def problem_2_10(self) -> List[Tuple[float, str, int]]:
        """Calculate the mass of a neutral atom of 235U in: a) amu, b) grams"""
        # a) Mass in amu
        m_235u = atomic_masses["235U"]
        # b) Mass in grams
        m_235u_g = m_235u / constants["N_A (1/(g-mol))"]
        return [(m_235u, "amu", 8), (m_235u_g, "g", 5)]

    def problem_2_11(self) -> str:
        """Show that 1 amu is numerically equal to the reciprocal of Avogadro's number"""
        return "Derivation Question"

    def problem_2_12(self) -> str:
        """Estimate the radius of the nucleus of 238U. Roughly what fraction of the 238U atom is
        taken up by the nucleus"""
        # The empirical atomic radius of uranium is 175 pm
        r = 175e-12  # m
        a = 238
        r_238u = self.formulae.nuclear_radius(a) / 1e15  # m
        # Fraction involves dividing the volume, simplified to:
        fraction = (r_238u / r) ** 3
        return f"{fraction * 100:.4g} %"

    def problem_2_13(self) -> List[Tuple[float, str, int]]:
        """Estimate the density of nuclear matter in g/cm^3 and kg/m^3. Take the mass of each
        nucleon to be 1.5*10^-24 g."""
        nucleon_mass = 1.5e-24
        # Suppose the number of nucleons is 10^39 and the nuclear radius is the coefficient of
        # the nuclear radius approximation formula:
        r, nucleons = 1.25, 1e39
        # The number of nucleons per cm^3 is:
        nucleon_density = 3 * nucleons / (4 * pi * r**3)
        # The density in g/cm^3:
        rho_g = nucleon_mass * nucleon_density
        # In kg/m^3:
        rho_kg = rho_g / 1000
        return [(rho_g, "g/cm^3", 3), (rho_kg, "kg/m^3", 3)]

    def problem_2_14(self) -> str:
        """The planet earth has a mass of approximately 6*10^24 kg. If the density of the earth
        were equal to to that of nuclei, how big would the earth be?"""
        approx_earth_mass = 6e24
        # Taking the result of the previous problem
        nuclear_density = self.problem_2_13()[1][0]
        # Solve for the radius in the relation rho = m / V
        radius = ((3 * approx_earth_mass) / (4 * pi * nuclear_density)) ** (1 / 3)
        return f"The earth's radius would be {radius:.3g} m if it was as dense as nuclear matter"

    def problem_2_15(self) -> Tuple[float, str, int]:
        """The complete combustion of 1 kg of bituminous coal releases about 3*10^7 J in heat
        energy. The conversion of 1 g of mass into energy is equivalent to the burning of how
        much coal?"""
        e = self.formulae.total_energy(1 / 1000)
        e_btu = e * conversion_factors["Btu/J"]
        # Hard coal is rated at 13000 Btu/lb and 1 kg = 2.205 lb. 1000 kg = 1 tonne
        mass_tonnes = e_btu / (13000 * 2.205 * 1000)
        return mass_tonnes, "tonnes", 3

    def problem_2_16(self) -> List[Tuple[float, str, int]]:
        """The fission of the nucleus of 235U releases approximated 200 MeV. How much energy in
        kW-hrs and MW-days is released when 1 g of 235U undergoes fission?"""
        a = 235
        n = constants["N_A (1/(g-mol))"] / a  # atoms/g
        e = n * 200  # MeV
        e_kwh = e * conversion_factors["J/MeV"] / conversion_factors["J/kWh"]
        # 1000 kW = 1 MW, 24 hrs = 1 day
        e_mwd = e_kwh / (1000 * 24)
        return [(e_kwh, "kWh", 3), (e_mwd, "MWd", 3)]

    def problem_2_17(self) -> Tuple[float, str, int]:
        """Compute the neutron-proton mass difference in MeV."""
        delta_m = abs(constants["m_p (MeV)"] - constants["m_n (MeV)"])
        return delta_m, "MeV", 4

    def problem_2_18(self) -> List[Tuple[float, str, int]]:
        """An electron starting from rest is accelerated across a potential difference of 5
        million volts.
        (a) What is its final kinetic energy?
        (b) What is its total energy?
        (c) What is its final mass?"""
        delta_e = 5e6
        q, m_e, c = constants["q_e (C)"], constants["m_e (kg)"], constants["c (m/s)"]
        # a) Since initial speed is 0, the kinetic energy is equal to the work done on the electron
        e_k = delta_e * q  # J
        # b) The total energy is the rest energy + the kinetic energy.
        e_rest = self.formulae.rest_energy(m_e)
        e_total = e_k + e_rest
        # c) e_total = m * c^2. The final mass calculation is simple
        m_f = e_total / c**2
        return [(e_k, "J", 3), (e_total, "J", 3), (m_f, "kg", 3)]

    def problem_2_19(self) -> str:
        """Derive Eq. (2.18)"""
        return "Derivation Question"

    def problem_2_20(self) -> str:
        """Show that the speed of any particle, relativistic or nonrelativistic, is given by the
        following formula: v = c * sqrt(1 - E_rest^2/E_total^2)"""
        return "Derivation Question"

    def problem_2_21(self) -> Tuple[float, str, int]:
        """Using the equation in Problem 2.20, calculated the speed of a 1 MeV electron, one with a
        kinetic energy of 1 MeV"""
        m_e, c = constants["m_e (kg)"], constants["c (m/s)"]
        e_rest = self.formulae.rest_energy(m_e) / conversion_factors["J/MeV"]
        e_total = e_rest + 1
        v = c * sqrt(1 - e_rest**2 / e_total**2)
        return v, "m/s", 3

    def problem_2_22(self) -> List[Tuple[float, str, int]]:
        """Compute the wavelengths of a 1 MeV: (a) photon, (b) neutron"""
        # a)
        e_t_photon = 1 * conversion_factors["J/MeV"]
        lambda_photon = self.formulae.photon_wavelength(e_t_photon)
        # b)
        m_n = constants["m_n (kg)"]
        lambda_neutron = self.formulae.wavelength_nonrel(m_n, e_t_photon)
        return [(lambda_photon, "m", 3), (lambda_neutron, "m", 3)]

    def problem_2_23(self) -> str:
        """Show that the wavelength of a relativistic particle is given by:
        lambda = lambda_C * m_e * c^2 / sqrt(e_total^2 - e_rest^2)"""
        return "Derivation Question"

    def problem_2_24(self) -> Tuple[float, str, int]:
        """Using the formula obtained in Problem 2.23, compute the wavelength of a 1 MeV electron"""
        m_e, h, c = constants["m_e (kg)"], constants["h (J-s)"], constants["c (m/s)"]
        e_rest = self.formulae.rest_energy(m_e)
        e_total = 1 * conversion_factors["J/MeV"] + e_rest
        lambda_c = h / (m_e * c)
        lambda_ = (lambda_c * m_e * c**2) / sqrt(e_total**2 - e_rest**2)
        return lambda_, "m", 3

    def problem_2_25(self) -> List[Tuple[float, str, int]]:
        """An electron moves with a kinetic energy equal to its rest-mass energy. Calculate the
        electron's:
        (a) total energy in units of m_e c^2
        (b) mass in units of m_e
        (c) speed in units of c
        (d) wavelength in units of the Compton wavelength"""
        m_e, h, c = constants["m_e (kg)"], constants["h (J-s)"], constants["c (m/s)"]
        # a)
        e_total = self.formulae.rest_energy(m_e) * 2
        # b)
        m = e_total / c**2
        # c)
        v = self.formulae.particle_speed(e_total, e_total / 2)
        # d)
        lambda_ = self.formulae.wavelength_w_compton(e_total, e_total / 2)
        lambda_c = h / (m_e * c)
        return [
            (e_total / (m_e * c**2), "m_e c^2", 3),
            (m / m_e, "m_e", 3),
            (v / c, "c", 3),
            (lambda_ / lambda_c, "lambda_C", 3),
        ]

    def problem_2_26(self) -> str:
        """A photon carries momentum, thus a free atom or nucleus recoils when it emits a photon.
        The energy of a photon is therefore less than the available transition energy (energy
        between states) by an amount equal to the recoil energy of the radiating system.
        (a) Given that E is the energy between two states and E_gamma is the energy of the emitted
        photon, show that E_gamma is approximately equal to E(1 - E / (2 * M * c^2)), where M is
        the mass of the atom or nucleus.
        (b) Compute E - E_gamma for the transitions from the first excited state of atomic hydrogen
        at 10.19 eV to ground and the first excited state of 12C at 4.43 MeV to ground
        """
        e_h, e_12c = (  # First excited states of H and 12C in J
            10.19e-6 * conversion_factors["J/MeV"],
            4.43 * conversion_factors["J/MeV"],
        )
        m_h, m_12c = (  # Masses of H and 12C in kg
            atomic_masses["H"]
            * conversion_factors["g/amu"]
            / conversion_factors["g/kg"],
            atomic_masses["12C"]
            * conversion_factors["g/amu"]
            / conversion_factors["g/kg"],
        )
        c = constants["c (m/s)"]
        e_gamma_h = e_h * (1 - e_h / (2 * m_h * c**2))
        e_gamma_12c = e_12c * (1 - e_12c / (2 * m_12c * c**2))
        delta_e_h = (e_h - e_gamma_h) / conversion_factors["J/MeV"] * 10e6
        delta_e_12c = (e_12c - e_gamma_12c) / conversion_factors["J/MeV"]
        return f"Derivation Question, {delta_e_h:.3g} eV, {delta_e_12c:.3g} MeV"

    def problem_2_27(self) -> str:
        """The first three excited states of the nucleus of 199Hg are at 0.158 MeV, 0.208 MeV, and
        0.403 MeV above the ground state. If all transitions between these states and ground state
        occured, what energy gamma-rays would be observed?"""
        e_0, e_1, e_2, e_3 = 0, 0.158, 0.208, 0.403
        # There are 6 possible transitions:
        gamma = [e_1 - e_0, e_2 - e_1, e_3 - e_2, e_3 - e_0, e_2 - e_0, e_3 - e_1]
        n = len(gamma)
        for i in range(n):
            gamma[i] = f"{gamma[i]:.3g} MeV"
        return ", ".join(gamma)

    def problem_2_28(self) -> str:
        """Using the chart of nuclides, complete the following reactions. If a daughter nucleus is
        also radioactive, indicate the complete decay chain.
        (a) 18N -> (b) 83Y -> (c) 135Sb -> (d) 219Rn ->"""
        result = []
        # a) 18N undergoes beta- decay. Therefore:
        result.append("18N -> 18O")
        # b) 83Y undergoes beta+ decay
        result.append("83Y -> 83Sr -> Rb -> 83Kr")
        # c) 135Sb undergoes beta- decay
        result.append("135Sb -> 135Te -> 135I -> 135Xe -> 135Cs -> 135Ba")
        # d) 219Rn undergoes alpha decay
        result.append("219Rn -> 215Po -> 211Pb -> 211Bi -> 207Tl -> 207Pb")
        return ", ".join(result)

    def problem_2_29(self) -> str:
        """Tritium (3H) decays by negative beta decay with a half-life of 12.26 years. The atomic
        weight of 3H is 3.016.
        (a) To what nucleus does H decay?
        (b) What is the mass in grams of 1mCi of tritium?"""
        # a) 3H undergoes B- decay to 3He, which is stable
        decay_chain = "3H -> 3He"
        # b)
        n_a, d = (
            constants["N_A (1/(g-mol))"],
            conversion_factors["Bq/Ci"],
        )
        m_3h, half_life = atomic_masses["3H"], half_lives["3H"]
        alpha = 1e-3 * d  # disintegrations/s
        lambda_ = self.formulae.decay_constant(half_life)
        num_atoms = alpha / lambda_
        mass = num_atoms * m_3h / n_a
        return f"{decay_chain}, {mass:.3g} g"

    def problem_2_30(self) -> Tuple[float, str, int]:
        """Approximately what mass of 90Sr (half-life = 28.8 yr) has the same activity as 1 g of
        60Co (half-life = 5.26 yr)?"""
        mass_60co = 1
        n_a = constants["N_A (1/(g-mol))"]
        t_sr, t_co = half_lives["90Sr"], half_lives["60Co"]
        m_90sr, m_60co = atomic_masses["90Sr"], atomic_masses["60Co"]
        # Compute the decay constants
        lambda_sr = self.formulae.decay_constant(t_sr)
        lambda_co = self.formulae.decay_constant(t_co)
        # Since the activities are equal, one can find the number of strontium atoms through the
        # equation: num_atoms_sr * lambda_sr = num_atoms_co * lambda_co
        # Where the number of atoms is given by m * N_A / M
        num_atoms_co = mass_60co * n_a / m_60co
        num_atoms_sr = num_atoms_co * lambda_co / lambda_sr
        mass_sr = num_atoms_sr * m_90sr / n_a
        return mass_sr, "g", 3

    def problem_2_31(self) -> Tuple[float, str, int]:
        """Carbon tetrachloride labeled with 14C is sold commercially with an activity of 10 mCi/mM.
        What fraction of the carbon atoms is 14C?"""
        n_a, d = (
            constants["N_A (1/(g-mol))"],
            conversion_factors["Bq/Ci"],
        )
        activity_ccl4, t_14c = 10e-3, half_lives["14C"]
        # CCl_4 has one carbon atom per molecule. Therefore:
        n_atoms_c = n_a * 1e-3  # atoms/mM
        activity_ccl4 *= d * 0.1  # disintegrations/s
        # Computing the decay constant and number of 14C atoms:
        lambda_ = self.formulae.decay_constant(t_14c)
        n_atoms_14c = activity_ccl4 / lambda_
        # The fraction is therefore:
        fraction = n_atoms_14c / n_atoms_c
        return fraction * 100, "%", 2

    def problem_2_32(self) -> str:
        """Titrated water (ordinary water containing some 1H 3HO) for some biological applications
        can be purchased in 1 cm^3 ampoules having an activity of 5 mCi/cm^3. What fraction of the
        water molecules contains an 3H atom?"""
        rho_1h3ho, m_1h3ho = 1.017, 20.04  # g/cm^3, g/mol
        n_a, d = (
            constants["N_A (1/(g-mol))"],
            conversion_factors["Bq/Ci"],
        )
        lambda_ = self.formulae.decay_constant(half_lives["3H"])
        v = m_1h3ho / rho_1h3ho  # cm^3/mol
        alpha = 5e-3 * v * d  # disintegrations/s
        n = alpha / lambda_
        fraction = (n / n_a) * 100
        return fraction, "%", 3

    def problem_2_33(self) -> List[Tuple[float, str, int]]:
        """After the initial cleanup effort at Three Mile Island, approximately 400,000 gallons of
        radioactive water remained in the basement of the containment building of the Three Mile
        Island Unit 2 nuclear plant. The principle sources of this radioactivity were 137Cs at 156
        uCi/cm^3 and 134Cs at 26 uCi/cm^3. How many atoms per cm^3 of these radionuclides were in
        the water at that time?"""
        alpha_137cs, alpha_134cs = 156, 26  # uCi/cm^3
        t_137cs, t_134cs = half_lives["137Cs"], half_lives["134Cs"]  # s
        d = conversion_factors["Bq/Ci"]
        # Decay constants
        lambda_137cs = self.formulae.decay_constant(t_137cs)
        lambda_134cs = self.formulae.decay_constant(t_134cs)
        # Number of atoms
        n_137cs = alpha_137cs * 1e-6 * d / lambda_137cs
        n_134cs = alpha_134cs * 1e-6 * d / lambda_134cs
        return [(n_137cs, "atoms/cm^3", 3), (n_134cs, "atoms/cm^3", 3)]

    def problem_2_34(self) -> List[Tuple[float, str, int]]:
        """One gram of 226Ra is placed in a sealed, evacuated capsule 1.2 cm^3 in volume.
        (a) At what rate does the helium pressure increase in the capsule, assuming all of the alpha
        particles are neutralized and retained in the free volume of the capsule?
        (b) What is the pressure 10 years after the capsule is sealed?"""
        n_a = constants["N_A (1/(g-mol))"]
        # a) The pressure rate is related to activity per unit volume:
        lambda_ = self.formulae.decay_constant(half_lives["226Ra"])
        n = 1 * n_a / (226 * 1.2)  # atoms/cm^3
        alpha = lambda_ * n  # disintegrations/(cm^3-s)
        # b) The activity per unit volume increases linearly with time
        t = 10 * 365 * 24 * 3600  # s
        p = alpha * t
        return [(alpha, "disintegrations/(cm^3-s)", 3), (p, "disintegrations/cm^3", 3)]

    def problem_2_35(self) -> Tuple[float, str, int]:
        """Polonium-210 decays to the ground state of 206Pb by the emission of a 5.305-MeV alpha
        particle with a half-life of 138 days. What mass of 21OPo is required to produce 1 MW of
        thermal energy from its radioactive decay?"""
        e, t = 5.305e6, 138 * 24 * 3600
        n_a, d = (
            constants["N_A (1/(g-mol))"],
            conversion_factors["Bq/Ci"],
        )
        # The power emitted per curie is:
        p = e * d * conversion_factors["J/MeV"] / conversion_factors["eV/MeV"]  # eV/Ci
        # To produce 1 MW, the desired activity is:
        alpha = 1e6 / p * d
        # Since alpha = lambda * N = lambda * m * N_A / M
        lambda_ = self.formulae.decay_constant(t)
        m = alpha * 210 / (lambda_ * n_a)
        return m / 1000, "kg", 4

    def problem_2_36(self) -> List[Tuple[float, str, int]]:
        """The radioisotope generator SNAP-9 was fueled with 475 g of 238PuC (plutonium-238
        carbide), which has a density of 12.5 g/cm^3. The 238 Pu has a half-life of 89 years and
        emits 5.6 MeV per disintegration, all of which may be assumed to be absorbed in the
        generator. The thermal to electrical efficiency of the system is 5.4 %. Calculate:
        (a) the fuel efficiency in curies per watt (thermal)
        (b) the specific power in watts (thermal) per gram of fuel
        (c) the power density in watts (thermal) per cm^3
        (d) the total electrical power of the generator"""
        m, rho, t, e = 475, 12.5, 89 * 365 * 24 * 3600, 5.6
        n_a = constants["N_A (1/(g-mol))"]
        # a)
        d = conversion_factors["Bq/Ci"]
        eff = 1 / e * 1 / d * 1 / conversion_factors["J/MeV"]
        # b)
        lambda_ = self.formulae.decay_constant(t)
        p = (
            lambda_
            * n_a
            / (
                rho
                * conversion_factors["L/mol"]
                * conversion_factors["mL/L"]
                * conversion_factors["Bq/Ci"]
                * eff
            )
        )
        # c)
        pd = rho * p
        # d)
        p_total = m * p * 5.4 / 100
        return [(eff, "Ci/W", 4), (p, "W/g", 4), (pd, "W/cm^3", 4), (p_total, "W", 4)]

    def problem_2_37(self) -> Tuple[float, str, int]:
        """Since the half-life of 235U (7.13 * 10^8 years) is less than that of 238U (4.51 * 10^9
        years), the isotopic abundance of 235U has been steadily decreasing since the earth was
        formed about 4.5 billion years ago. How long ago was the isotopic abundance of 235U equal
        to 3.0 a/o, the enrichment of the uranium used in many nuclear power plants?"""
        # Compute decay constants
        lambda_238 = self.formulae.decay_constant(half_lives["238U"])
        lambda_235 = self.formulae.decay_constant(half_lives["235U"])
        # Abundances
        x_235 = 0.03
        x_238 = 1 - x_235
        # The primary abundance of 235U is 0.72 a/o. Computing the time t, after rearranging the
        # decay equation:
        a = 0.72 / 100
        t = abs(log(a * x_238 / (x_235 * (1 - a))) / (lambda_235 - lambda_238))
        t /= 365 * 24 * 3600
        return t, "yr", 3

    def problem_2_38(self) -> str:
        """The radioactive isotope Y is produced at the rate of R atoms/sec by neutron bombardment
        of X according to the reaction X(n, gamma)Y. If the neutron bombardment is carried out for a
        time equal to the half-life of Y, what fraction of the saturation activity of Y will be
        obtained assuming that there is no Y present at the start of the bombardment?"""
        return "Derivation Question"

    def problem_2_39(self) -> str:
        """Consider the chain decay with no atoms of B present at t = O:
        A -> B -> C ->
        (a) Show that the activity of B rises to a maximum value at the time t_m given by:
        t_m = 1/(lambda_B - lambda_A) * ln(lambda_B/lambda_A)
        at which time the activities of A and B are equal.
        (b) Show that, for t  < t_m, the activity of B is less than that of A, whereas the reverse
        is the case for t > t_m."""
        return "Derivation Question"

    def problem_2_40(self) -> str:
        """Show that if the half-life of B is much shorter than the half-life of A, then the
        activities of A and B in Problem 2.39 eventually approach the same value. In this case, A
        and B are said to be in secular equilibrium."""
        return "Derivation Question"

    def problem_2_41(self) -> str:
        """Show that the abundance of 234U can be explained by assuming that this isotope originates
        solely from the decay of 238U"""
        return "The half-life of 238U >> others in the decay chain, it has the most influence."

    def problem_2_42(self) -> Tuple[float, str, int]:
        """Radon-222, a highly radioactive gas with a half-life of 3.8 days that originates in the
        decay of 234U (see the chart of nuclides), may be present in uranium mines in dangerous
        concentrations if the mines are not properly ventilated. Calculate the activity of 222Rn in
        Bq per metric ton of natural uranium"""
        gamma_i = [abundances["234U"], abundances["235U"], abundances["238U"]]
        m_i = [atomic_masses["234U"], atomic_masses["235U"], atomic_masses["238U"]]
        # Since the half-life of 234U is much greater than that of half-life 222Rn, the activities
        # are approximately equal. Computing the abundance of 234U
        gamma = self.formulae.abundance_natural(gamma_i, m_i, gamma_i[0], m_i[0])
        # Obtaining the activity in disintegrations/s
        lambda_ = self.formulae.decay_constant(half_lives["234U"])
        alpha = self.formulae.activity_atom_density(lambda_, m_i[0], gamma)
        alpha *= 1e6 / conversion_factors["Bq/Ci"]
        return alpha, "Ci/ton", 3

    def problem_2_43(self) -> str:
        """According to U.S. Nuclear Regulatory Commission regulations, the maximum permissible
        concentration of radon-222 in air in equilibrium with its short-lived daughters is 3 pCi/L
        for nonoccupational exposure. This corresponds to how many atoms of radon-222 per cm^3?
        """
        # Since 1 mL = 1 cm^3, the activity can be represented as:
        alpha = (
            3e-12 / conversion_factors["mL/L"] * conversion_factors["Bq/Ci"]
        )  # (disintegrations/s)/cm^3
        # The short-lived daughters are 218Po, 214Pb, 214Bi, and 214Po based on the decay chain.
        t_218po, t_214po = half_lives["218Po"], half_lives["214Po"]
        t_214pb, t_214bi = half_lives["214Pb"], half_lives["214Bi"]
        # Computing the total half-life of short-lived daughters:
        t_total = sum([t_214po, t_218po, t_214bi, t_214pb])
        # Since alpha = lambda * N, N = alpha/lambda
        lambda_ = self.formulae.decay_constant(t_total)
        n = alpha / lambda_
        return n, "atoms/cm^3", 3

    def problem_2_44(self) -> str:
        """Consider again the decay chain in Problem 2.39 in which the nuclide A is produced at the
        constant rate of R atoms/s. Derive an expression for the activity of B as a function of
        time."""
        return "Derivation Question"

    def problem_2_45(self) -> List[Tuple[float, str, int]]:
        """Complete the following reactions and calculate their Q values. [Note: The atomic weight
        of 14C is 14.003242.]
        (a) 4He(p,d) (b) 9Be(alpha,n) (c) 14N(n,p) (d) 115In(d,p) (e) 207Pb(gamma,n)"""
        # a) The reaction is 4He(p,d)3He or 4He + 1H -> 3He + 2H
        q_a = self.formulae.q_value_atomic_masses(
            atomic_masses["4He"],
            atomic_masses["1H"],
            atomic_masses["3He"],
            atomic_masses["2H"],
        )
        # b) The reaction is 9Be(alpha,n)12C or 9Be + 4He -> 12C + n
        q_b = self.formulae.q_value_atomic_masses(
            atomic_masses["9Be"],
            atomic_masses["4He"],
            atomic_masses["12C"],
            atomic_masses["n"],
        )
        # c) The reaction is 14N(n,p)14C or 14N + n -> 14C + 1H
        q_c = self.formulae.q_value_atomic_masses(
            atomic_masses["14N"],
            atomic_masses["n"],
            atomic_masses["14C"],
            atomic_masses["1H"],
        )
        # d) The reaction is 115In(d,p)116In or 115In + 2H -> 116In + 1H
        q_d = self.formulae.q_value_atomic_masses(
            atomic_masses["115In"],
            atomic_masses["2H"],
            atomic_masses["116In"],
            atomic_masses["1H"],
        )
        # e) The reaction is 207Pb(gamma, n)206Pb or 207Pb + gamma -> 206Pb + n
        q_e = self.formulae.q_value_atomic_masses(
            atomic_masses["207Pb"],
            atomic_masses["gamma"],
            atomic_masses["206Pb"],
            atomic_masses["n"],
        )
        return [(i, "MeV", 5) for i in [q_a, q_b, q_c, q_d, q_e]]

    def problem_2_46(self) -> List[Tuple[float, str, int]]:
        """
        (a) Compute the recoil energy of the residual, daughter nucleus following the emission of a
        4.782-MeV alpha-particle by 226Ra.
        (b) What is the total disintegration energy for this decay process?
        """
        e_alpha = 4.782  # MeV
        # a) The decay chain is 226Ra -> 222Rn + alpha (4He). The recoil energy is the difference
        # between the Q value of the reaction and emission energy
        q = self.formulae.q_value_atomic_masses(
            atomic_masses["226Ra"],
            0,
            atomic_masses["222Ra"],
            atomic_masses["4He"],
        )
        recoil = q - e_alpha
        # b) The disintegration energy is equal to the Q value
        return [(i, "MeV", 4) for i in [recoil, q]]

    def problem_2_47(self) -> str:
        """In some tabulations, atomic masses are given in terms of the mass excess rather than as
        atomic masses. The mass excess, Δ, is the difference Δ = M - A, Where M is the atomic mass
        and A is the atomic mass number. For convenience, Δ, which may be positive or negative, is
        usually given in units of MeV. Show that the Q value is given by Q = (Δa + Δb) - (Δc + Δd)
        """
        # See formula q_value_mass_defects() for the result of this derivation and the next problem
        return "Derivation Question"

    def problem_2_48(self) -> Tuple[float, str, int]:
        """The mass excesses for the (neutral) atoms in the reaction in Example 2.8 are as follows:
        Δ(3H) = 14.95 MeV, Δ(2H) = 13.14 MeV, Δ(n) = 8.07 MeV, and Δ(4He) = 2.42 MeV. Calculate the
        Q value of this reaction using the results of Problem 2.47"""
        # The reaction is 3H + 2H -> 4He + n. Therefore:
        q = self.formulae.q_value_mass_defects(14.95, 13.14, 8.07, 2.42)
        return q, "MeV", 5

    def problem_2_49(self) -> Tuple[float, str, int]:
        """The atomic weight of 206Pb is 205.9745. Using the data in Problem 2.35, calculate the
        atomic weight of 210Po"""
        e_alpha = 5.305  # The alpha particle energy from Problem 2.35
        # The reaction is 210Po -> 206Pb + alpha. The recoil energy is assumed to be 0 because of
        # how much heavier the nuclide is compared to the alpha particle. Therefore Q = E_alpha.
        q = e_alpha
        # Rearranging the Q value formula:
        m_210po = (
            q / (conversion_factors["MeV/amu"])
            + atomic_masses["206Pb"]
            + atomic_masses["4He"]
        )
        return m_210po, "amu", 7

    def problem_2_50(self) -> List[Tuple[float, str, int]]:
        """Tritium (3H) can be produced through the absorption of low-energy neutrons by deutrerium.
        The reaction is 2H + n → 3H + gamma, Where the gamma-ray has an energy of 6.256 MeV.
        (a) Show that the recoil energy of the 3H nucleus is approximately 7 KeV.
        (b) What is the Q-value of the reaction?
        (c) Calculate the separation energy of the of the last neutron in 3H
        (d) Using the binding energy of 2H of 2.23 MeV and the result from part (c), compute the
        total binding energy of 3H
        """
        e_gamma, be_2h = 6.256, 2.23  # MeV
        m_2h, m_3h, m_n = atomic_masses["2H"], atomic_masses["3H"], atomic_masses["n"]
        # a) The recoil energy can be obtained from the relation between kinetic energy and
        # momentum: 1/2 mv^2 = p^2/(2M), since p = mv = E/c: K = E^2/(2MC^2)
        recoil = e_gamma**2 / (2 * m_3h * conversion_factors["MeV/amu"])
        recoil *= 1000  # keV
        # b) Q-value
        q = self.formulae.q_value_atomic_masses(m_2h, m_n, 0, m_3h)
        # c) E_s = binding energy of the last neutron
        e_s = self.formulae.separation_energy(m_2h, m_3h)
        # d) The total binding energy is the separation energy of the last neutron plus the binding
        # energy of the ^{A-1}H nuclide, which is BE(2H)
        be = e_s + be_2h
        return [(recoil, "keV", 4), (q, "MeV", 5), (e_s, "MeV", 5), (be, "MeV", 5)]

    def problem_2_51(self) -> str:
        """Consider the reaction 6Li(alpha,p)9Be. Using atomic mass data, compute:
        (a) the total binding energy of 6Li, 9Be, and 4He;
        (b) the Q-value of the reaction using the results in part (a)"""
        m_6li, m_9be, m_4he = (
            atomic_masses["6Li"],
            atomic_masses["9Be"],
            atomic_masses["4He"],
        )
        # a) The binding energies are the separation energies of the last neutrons
        be_6li = self.formulae.binding_energy(m_6li, 3, 3)
        be_9be = self.formulae.binding_energy(m_9be, 4, 5)
        be_4he = self.formulae.binding_energy(m_4he, 2, 2)
        # b) The Q-value can be obtained with atomic masses. Note: Alpha particle mass is
        # approximately equal to M(4H)
        q = self.formulae.q_value_atomic_masses(
            m_6li, atomic_masses["alpha"], atomic_masses["p"], m_9be
        )
        return ", ".join(f"{i:.5g} MeV" for i in [be_6li, be_9be, be_4he, q])

    def problem_2_52(self) -> List[Tuple[float, str, int]]:
        """Using atomic mass data, compute the average binding energy per nucleon of the following
        nuclei: (a) 2H (b) 4He (c) 12C (d) 51V (e) 138Ba (f) 235U"""
        # Divide the binding energies by the atomic number
        be_2h = self.formulae.binding_energy(atomic_masses["2H"], z=1, n=1) / 2
        be_4he = self.formulae.binding_energy(atomic_masses["4He"], z=2, n=2) / 4
        be_12c = self.formulae.binding_energy(atomic_masses["12C"], z=6, n=6) / 12
        be_51v = self.formulae.binding_energy(atomic_masses["51V"], z=23, n=28) / 51
        be_138ba = (
            self.formulae.binding_energy(atomic_masses["138Ba"], z=56, n=82) / 138
        )
        be_235u = self.formulae.binding_energy(atomic_masses["235U"], z=92, n=143) / 235
        return [
            (i, "MeV", 4) for i in [be_2h, be_4he, be_12c, be_51v, be_138ba, be_235u]
        ]

    def problem_2_53(self) -> List[Tuple[float, str, int]]:
        """Using the mass formula, compute the binding energy per nucleon for the nuclei in Problem
        2.52. Compare the results with those obtained in that problem"""
        be_2h = self.formulae.binding_energy_mass_eqn(n=1, a=2, z=1) / 2
        be_4he = self.formulae.binding_energy_mass_eqn(n=2, a=4, z=2) / 4
        be_12c = self.formulae.binding_energy_mass_eqn(n=6, a=12, z=6) / 12
        be_51v = self.formulae.binding_energy_mass_eqn(n=28, a=51, z=23) / 51
        be_138ba = self.formulae.binding_energy_mass_eqn(n=82, a=138, z=56) / 138
        be_235u = self.formulae.binding_energy_mass_eqn(n=143, a=235, z=92) / 235
        return [
            (i, "MeV", 4) for i in [be_2h, be_4he, be_12c, be_51v, be_138ba, be_235u]
        ]

    def problem_2_54(self) -> List[Tuple[float, str, int]]:
        """Compute the separation energies of the last neutron in the following nuclei:
        (a) 4He (b) 7Li (c) 17O (d) 51V (e) 208Pb (f) 235U"""
        # Using the separation energy formula:
        e_4he = self.formulae.separation_energy(
            atomic_masses["3He"], atomic_masses["4He"]
        )
        e_7li = self.formulae.separation_energy(
            atomic_masses["6Li"], atomic_masses["7Li"]
        )
        e_17o = self.formulae.separation_energy(
            atomic_masses["16O"], atomic_masses["17O"]
        )
        e_51v = self.formulae.separation_energy(
            atomic_masses["50V"], atomic_masses["51V"]
        )
        e_208pb = self.formulae.separation_energy(
            atomic_masses["207Pb"], atomic_masses["208Pb"]
        )
        e_235u = self.formulae.separation_energy(
            atomic_masses["234U"], atomic_masses["235U"]
        )
        return [(i, "MeV", 4) for i in [e_4he, e_7li, e_17o, e_51v, e_208pb, e_235u]]

    def problem_2_55(self) -> str:
        """Derive Eq.(2.53)"""
        return "Derivation Question"

    def problem_2_56(self) -> Tuple[float, str, int]:
        """What is 1 atmosphere pressure in units of eV/cm^3?"""
        # Since 1 Pa = 1 N/m^2 and 1 J = 1 N-m and the conversion of MeV to eV and m to cm^3 cancels
        # out in terms of numerical value:
        p = 1 * conversion_factors["Pa/atm"] / conversion_factors["J/MeV"]
        return p, "eV/cm^3", 4

    def problem_2_57(self) -> Tuple[float, str, int]:
        """Calculate the atom density of graphite having a density of 1.60 g/cm^3"""
        # Graphite is comprised of carbon, therefore:
        n = self.formulae.atom_density(rho=1.6, m=atomic_masses["C"])
        return n, "atoms/cm^3", 3

    def problem_2_58(self) -> Tuple[float, str, int]:
        """Calculate the activity of 1 gram of natural uranium"""
        # The average atomic weight of naturally occuring uranium was computed in Problem 2.6:
        m = self.problem_2_6()[0]
        # The activity is therefore the sum of the activities of 234U, 235U, and 238U:
        t = [half_lives["234U"], half_lives["235U"], half_lives["238U"]]
        gamma = [abundances["234U"], abundances["235U"], abundances["238U"]]
        n, alpha = len(t), 0
        for i in range(n):
            lambda_ = self.formulae.decay_constant(t[i])
            alpha += self.formulae.activity_atom_density(lambda_, m, gamma[i])
        alpha /= conversion_factors["Bq/Ci"]
        return alpha, "Ci", 4

    def problem_2_59(self) -> Tuple[float, str, int]:
        """What is the atom density of 235U in uranium enriched to 2.5 a/o in this isotope if the
        physical density of the uranium is 19.0 g/cm^3"""
        # The average atomic weight of naturally occuring uranium was computed in Problem 2.6:
        m, rho, w = self.problem_2_6()[0], 19.0, 2.5
        n = self.formulae.atom_density_component(w, rho, m)
        return n, "atoms/cm^3", 3

    def problem_2_60(self) -> List[Tuple[float, str, int]]:
        """Plutonium-239 undergoes alpha-decay with a half-life of 24,000 years. Compute the
        activity of 1 gram of plutonium dioxide, 239PuO_2. Express the activity in terms of Ci
        and Bq"""
        m_239pu = atomic_masses["239Pu"]
        m_total = m_239pu + 2 * atomic_masses["O"]
        gamma = m_239pu / m_total
        lambda_ = self.formulae.decay_constant(half_lives["239Pu"])
        alpha_bq = self.formulae.activity_atom_density(lambda_, m_239pu, gamma)
        alpha_ci = alpha_bq / conversion_factors["Bq/Ci"]
        return [(alpha_ci, "Ci", 4), (alpha_bq, "Bq", 4)]

    def problem_2_61(self) -> List[Tuple[float, str, int]]:
        """It has been proposed to use uranium carbide (UC) for the initial fuel in certain types of
        breeder reactors, with the uranium enriched to 25 w/o. The density of UC is 13.6 g/cm^3
        (a) What is the atomic weight of the uranium?
        (b) What is the atom density of the 235U?"""
        # a) Found by computing the average atomic weight
        m_u = self.formulae.average_atomic_weight(
            gamma=[25, 75], m=[atomic_masses["235U"], atomic_masses["238U"]]
        )
        # b) Computing the weight percent of uranium, then the atom density of uranium-235 in UC
        w = self.formulae.weight_percent(x=m_u, m=1, y=atomic_masses["C"], n=1)
        rho_u = 13.6 * w / 100
        n = self.formulae.atom_density_component(
            w=25, rho=rho_u, m=atomic_masses["235U"]
        )
        return [(m_u, "amu", 4), (n, "atoms/cm^3", 4)]

    def problem_2_62(self) -> List[Tuple[float, str, int]]:
        """Compute the atom densities of 235U and 238U in UO_2 of physical density 10.8 g/cm^3 if
        the uranium is enriched to 3.5 w/o in 235U"""
        m_235u, m_238u = atomic_masses["235U"], atomic_masses["238U"]
        m_u = self.formulae.average_atomic_weight(gamma=[3.5, 96.5], m=[m_235u, m_238u])
        m_uo_2 = m_u + 2 * atomic_masses["O"]
        w = m_u / m_uo_2
        rho_u = 10.8 * w
        n_235u = self.formulae.atom_density_component(w=3.5, rho=rho_u, m=m_235u)
        n_238u = self.formulae.atom_density_component(w=96.5, rho=rho_u, m=m_238u)
        return [(n_235u, "atoms/cm^3", 4), (n_238u, "atoms/cm^3", 4)]

    def problem_2_63(self) -> List[Tuple[float, str, int]]:
        """The fuel for a certain breeder reactor consists of pellets composed of mixed oxides, UO_2
        and PuO_2, with the PuO_2 comprising approximately 30 w/o of the mixture. The uranium is
        essentially all 238U, whereas the plutonium contains the following isotopes: 239Pu
        (70.5 w/o), 240Pu(21.3 w/o), 241Pu(5.5 w/o), and 242Pu(2.7 w/o). Calculate the number of
        atoms of each isotope per gram of the fuel."""
        # Computing the average atomic weight
        gamma = [70.5, 21.3, 5.5, 2.7]
        m = [
            atomic_masses["239Pu"],
            atomic_masses["240Pu"],
            atomic_masses["241Pu"],
            atomic_masses["242Pu"],
        ]
        m_avg = self.formulae.average_atomic_weight(gamma, m)
        m_pu_o_2 = m_avg + 2 * atomic_masses["O"]
        w = m_avg / m_pu_o_2
        rho = 0.3
        j, n = len(gamma), []
        # Iteratively calculating the atom densities for each isotope
        for i in range(j):
            n.append(gamma[i] * self.formulae.atom_density_component(w, rho, m[i]))
        return [(k, "atoms/g", 4) for k in n]


if __name__ == "__main__":
    ch2p = Chapter2Problems()
    ans = AnswerWriter(obj=ch2p, num_questions=63, path="./answers/problems.txt")
    ans.write_answers()
