"""
chapter_3_examples.py

"""

import math
from helpers.answer_types import SingleAnswer, MultiPartAnswer
from helpers import conversions, data
from helpers.formulae import Formulae
from helpers.answer_writer import AnswerWriter

formulae = Formulae()


class Chapter3Examples:
    """Solutions of Chapter 3 example problems"""

    type_and_chapter = ("Example", 3)

    def example_3_1(self) -> MultiPartAnswer:
        """A beam of 1-MeV neutrons of intensity 5 * 10^8 neutrons/(cm^2-s) strikes a thin 12C
        target. The area of the target is 0.5 cm^2 and is 0.05 cm thick. The beam has a
        cross-sectional area of 0.1 cm^2 At 1 MeV, the total cross-section of 12C is 2.6 b.
        (a) At what rate do interactions take place in the target?
        (b) What is the probability that a neutron in the beam will have a collision in the target?
        """
        sigma_t = 2.6 / conversions.area["b/cm^2"]
        collisions = formulae.ch3.se_energy_beam_collisions(
            sigma=sigma_t, i=5e8, n=data.atom_densities["C"], a=0.1, x=0.05
        )  # interactions/s
        neutron_strikes = 5e8 * 0.1  # neutrons/s
        # Divide the number of collisions by the number of neutron strikes
        p = collisions / neutron_strikes * 100
        return MultiPartAnswer(
            ans=[collisions, p], units=["interactions/s", "%"], sig_figs=3
        ).format()

    def example_3_2(self) -> SingleAnswer:
        """There are only two absorption reactions - namely, radiative capture and fission - that
        can occur when 0.0253-eV neutrons interact with 235U. The cross-sections for these reactions
        are 99 b and 582 b, respectively. When a 0.0253-e V neutron is absorbed by 235U, what is the
        relative probability that fission will occur?"""
        # Since the radiative capture cross-section and fission cross-section are proportional to
        # the probabilities of radiative capture and fission:
        sigma_gamma, sigma_f = 99, 582
        p = sigma_f / (sigma_gamma + sigma_f) * 100
        return SingleAnswer(ans=p, units="%", sig_figs=3).format()

    def example_3_3(self) -> MultiPartAnswer:
        """Referring to Example 3.1, calculate the
        (a) macroscopic total cross-section of 12C at 1-MeV;
        (b) collision density in the target."""
        sigma_t = 2.6 / conversions.area["b/cm^2"]
        n, i = data.atom_densities["C"], 5e8
        sigma_macro_t = formulae.ch3.macroscopic_cross_section(n, sigma_t)
        f = formulae.ch3.collision_density(i, n, sigma_t)
        return MultiPartAnswer(
            ans=[sigma_macro_t, f], units=["cm^-1", "collisions/(cm^3-s)"], sig_figs=3
        ).format()

    def example_3_4(self) -> SingleAnswer:
        """Calculate the mean free path of 100-keV neutrons in liquid sodium. At this energy, the
        total cross-section of sodium is 3.4 b."""
        sigma_t, n = 3.4 / conversions.area["b/cm^2"], data.atom_densities["Na"]
        sigma_macro_t = formulae.ch3.macroscopic_cross_section(n, sigma_t)
        lambda_ = formulae.ch3.mean_free_path(sigma_macro_t)
        return SingleAnswer(ans=lambda_, units="cm", sig_figs=3).format()

    def example_3_5(self) -> SingleAnswer:
        """The absorption cross-section of 235U and 238U at 0.0253 eV are 680.8 b and 2.70 b,
        respectively. Calculate the macroscopic cross-section for natural uranium at this energy.
        """
        sigma_235u, sigma_238u = (
            680.8 / conversions.area["b/cm^2"],
            2.7 / conversions.area["b/cm^2"],
        )
        # From Chapter 2, the atom densities are:
        n_235u, n_238u = 3.48e20, 4.83e22
        sigma_macro = formulae.ch3.macroscopic_cross_section_mixture(
            n_235u, sigma_235u, n_238u, sigma_238u
        )
        return SingleAnswer(ans=sigma_macro, units="cm^-1", sig_figs=3).format()

    def example_3_6(self) -> MultiPartAnswer:
        """The scattering cross-sections (in barns) of hydrogen and oxygen at 1 MeV and 0.0253 eV
        are given in the following table.
           | 1 MeV | 0.0253 eV |
        ---+-------+-----------+
         H |   3   |    21     |
         O |   8   |     4     |
        What are the values of as for the water molecule at these energies?"""
        # At 1 MeV, one gets:
        macro_sigma_1 = formulae.ch3.molecule_cross_section(
            m=2, sigma_m=3, n=1, sigma_n=8
        )
        # This equation does not hold at 0.0253 eV, and the experimental value is 103
        macro_sigma_2 = 103
        return MultiPartAnswer(
            ans=[float(macro_sigma_1), float(macro_sigma_2)], units="b", sig_figs=3
        ).format()

    def example_3_7(self) -> SingleAnswer:
        """A certain research reactor has a flux of 1 * 10^13 neutrons/(cm^2-s) and a volume of
        64,000 cm^3. If the macroscopic fission cross-section in the reactor is 0.1 cm^-1, what is
        the power of the reactor?"""
        # The energy released in fission is 200 MeV
        phi, volume, sigma_macro_f, energy = 1e13, 64e3, 0.1, 200
        # Computing the power rate
        power_rate = energy * conversions.energy["J/MeV"] * 1e-6  # MW/(fission/s)
        # The fission rate is:
        f = formulae.ch3.collision_density_flux(sigma_macro_f, phi)
        # Power density
        pd = power_rate * f  # MW/cm^3
        # The total power is the power density multiplied by the reactor volume
        p = pd * volume
        return SingleAnswer(ans=p, units="MW", sig_figs=3).format()

    def example_3_8(self) -> SingleAnswer:
        """Using experimental elastic scattering data, estimate the radius of the C nucleus."""
        # From Fig. 3.4, the scattering is 4.8 b from 0.02 eV to 0.01 MeV due to potential
        # scattering. Therefore:
        sigma_e = 4.8 / conversions.area["b/cm^2"]
        r = math.sqrt(sigma_e / (4 * math.pi))
        return SingleAnswer(ans=r, units="cm", sig_figs=3).format()

    def example_3_9(self) -> SingleAnswer:
        """The value of the radiative capture cross-section for 1H at 0.0253 eV is 0.332 b. What is
        the radiative capture cross-section at 1 eV?"""
        # Since sigma_gamma is 1/v, it can be written as
        # sigma_gamma(E) = sigma_gamma(E_0) * sqrt(E/E_0)
        # where E_0 is any energy. Therefore:
        sigma_gamma_1 = 0.332
        sigma_gamma_2 = sigma_gamma_1 * math.sqrt(0.0253 / 1)
        return SingleAnswer(ans=sigma_gamma_2, units="b", sig_figs=3).format()

    def example_3_10(self) -> MultiPartAnswer:
        """A 1 MeV neutron is scattered through an angle of 45 deg in a collision with a 2H nucleus.
        (a) What is the energy of the scattered neutron?
        (b) What is the energy of the recoiling nucleus?
        (c) How much of a change in lethargy does the neutron undergo in this collision?
        """
        theta, a, e_i = 45 / conversions.angle["deg/rad"], 2, 1
        # a)
        e_s = formulae.ch3.energy_after_collision(e_i, a, theta)
        # b) Recoil energy is the initial energy minus the scattering energy
        e_r = e_i - e_s
        # c)
        delta_u = abs(formulae.ch3.lethargy(e_i, e_i) - formulae.ch3.lethargy(e_i, e_s))
        return MultiPartAnswer(
            ans=[e_s, e_r, delta_u], units=["MeV", "MeV", ""], sig_figs=3
        ).format()

    def example_3_11(self) -> SingleAnswer:
        """A small indium foil is placed at a point in a reactor where the 2,200 meters-per second
        flux is 5 * 10^12 neutrons/(cm^2-s). The neutron density can be represented by a Maxwellian
        function with a temperature of 600 deg C. At what rate are the neutrons absorbed per cm^3 in
        the foil?"""
        # Indium is a non-1/v absorber, therefore:
        sigma_macro_a, phi = data.cross_sections["In"]["macro_a"], 5e12
        g = data.non_1_v_factors["In"]["g_a"]["600C"]
        f = formulae.ch3.absorption_rate_non_1_v_alt(g, sigma_macro_a, phi)
        return SingleAnswer(ans=f, units="neutrons/cm^3", sig_figs=3).format()

    def example_3_12(self) -> SingleAnswer:
        """The total initial fuel loading of a particular reactor consists of 120 fuel rods. After
        the reactor has been operated at a steady power of 100 MW for 1 year, the fuel is removed.
        Assuming that all rods contribute equally to the total power, estimate the activity of a
        fuel rod 1 day after removal."""
        # Converting the operating time to days:
        operating_time = 1 * conversions.time["s/yr"] / conversions.time["s/day"]
        p, removal_time = 100, 1
        # Dividing by 120 rods yields the answer:
        alpha = formulae.ch3.fission_product_activity_reactor(
            p, operating_time, removal_time
        )
        # Dividing by 120 rods yields the answer:
        alpha /= 120
        return SingleAnswer(ans=alpha, units="Ci", sig_figs=2).format()

    def example_3_13(self) -> SingleAnswer:
        """Calculate the value of eta for natural uranium at 0.0253 eV."""
        # 238U does not fission with low-energy neutrons. Therefore the numerator of the equation
        # used below only involves values for 235U. Since the atom density is proportional to the
        # isotopic abundance:
        nu_235u, sigma_f_235u, sigma_a_235u, gamma_235u = (
            data.fissile_nuclide_thermal_data["235U"]["nu"],
            data.fissile_nuclide_thermal_data["235U"]["f"],
            data.fissile_nuclide_thermal_data["235U"]["a"],
            data.abundances["235U"],
        )
        sigma_a_238u, gamma_238u = (
            data.nuclide_cross_sections["238U"]["a"],
            data.abundances["238U"],
        )
        eta = (
            nu_235u
            * gamma_235u
            * sigma_f_235u
            / (gamma_235u * sigma_a_235u + gamma_238u * sigma_a_238u)
        )
        # Note: For reasons discussed in Chap. 6, this is not precisely the value used in reactor
        # problems involving natural uranium; see Example 6.11.
        return SingleAnswer(ans=eta, units="", sig_figs=4).format()

    def example_3_14(self) -> MultiPartAnswer:
        """The energy released by the fissioning of 1 g of 235U is equivalent to the combustion of
        how much
        (a) coal with a heat content of 3 x 10^7 J/kg (13,000 Btu/lb)
        (b) oil at 4.3 x 10^7 J/kg (6.5 x 106 Btu/barrel)?"""
        # The fissioning of 1 g of uranium-235 is equal to approximately 1 MW/day = 24000 kWh
        # = 8.64 x 10^10 J.
        p = 8.64e10  # J
        # a)
        coal_heat = 3e7  # J/kg
        m_1 = p / coal_heat
        # b)
        oil_heat = 4.3e7  # J/kg
        m_2 = p / oil_heat
        formulae.ch3.scattered_photon_energy(e=1, theta=1)
        return MultiPartAnswer(ans=[m_1, m_2], units="kg", sig_figs=3).format()

    def example_3_15(self) -> MultiPartAnswer:
        """Calculate the mass attenuation coefficient of UO_2 for 1 MeV gamma-rays. What is their
        mean free path? The density of UO_2 is about 10 g/cm^3"""
        m_u, m_o, rho = data.atomic_masses["U"], data.atomic_masses["O"], 10
        # Calculating the molecular mass and weight percents
        m = m_u + 2 * m_o
        omega_u = m_u / m * 100
        omega_o = 100 - omega_u
        # From Table II.4:
        mu_star_u = data.mass_attenuation_coefficients["U"]["1.0 MeV"]
        mu_star_o = data.mass_attenuation_coefficients["O"]["1.0 MeV"]
        # Calculating the mass attenuation coefficient of the mixture:
        mu_star = formulae.ch3.mass_attenuation_coeff_mixture(
            omega=[omega_u, omega_o], mu_star=[mu_star_u, mu_star_o]
        )
        # The mean free path is:
        mu = mu_star * rho
        lambda_ = formulae.ch3.gamma_ray_mean_free_path(mu)
        return MultiPartAnswer(
            ans=[mu_star, lambda_], units=["g/cm^2", "cm"], sig_figs=3
        ).format()

    def example_3_16(self) -> SingleAnswer:
        """It is proposed to store liquid radioactive waste in a steel container. If the intensity
        of gamma-rays incident on the interior surface of the tank is estimated to be 3 * 10^11
        gamma-rays/(cm^2-s) and the average gamma-ray energy is 0.8 MeV, at what rate is energy
        deposited at the surface of the container?"""
        i, e = 3e11, 0.8  # gamma-rays/(cm^2-s), MeV
        # Steel is mostly iron. Therefore:
        mu_rho_a = data.mass_absorption_coefficients["Fe"]["0.8 MeV"]
        w = formulae.ch3.energy_deposition_rate_att(e, i, mu_rho_a)
        return SingleAnswer(ans=w, units="MeV/(g-s)", sig_figs=3).format()

    def example_3_17(self) -> SingleAnswer:
        """Sodium-24 (T_1/2 = 15 hr) is often used in medicine as a radioactive tracer. It emits
        beta-rays with a maximum energy of 1.39 MeV. What is the maximum range of the beta-rays in
        animal tissue?"""
        e_max = 1.39
        # Since the density of most animal tissue approximately 1 g/cm:
        r = formulae.ch3.beta_ray_max_range(e_max, rho=1)
        return SingleAnswer(ans=r, units="cm", sig_figs=3).format()


if __name__ == "__main__":
    ch3e = Chapter3Examples()
    writer = AnswerWriter(obj=ch3e, num_questions=17, path="./answers/examples.txt")
    writer.write_answers()
