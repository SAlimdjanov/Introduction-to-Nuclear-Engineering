"""
chapter_3.py

"""

import math
from typing import List, Union
from .. import conversions, data


class RadiationMatterInteraction:
    """Chapter 3"""

    def se_beam_intensity(self, n: float, v: float) -> float:
        """Computes the intensity of a single-energy neutron beam given the density and speed of
        the neutrons"""
        return n * v

    def se_energy_beam_collisions(
        self, sigma: float, i: float, n: float, a: float, x: float
    ) -> float:
        """Calculates the number of collisions per second in a single-energy neutron beam given the
        proportionality constant, intensity, target atom density, target area, and target thickness
        """
        return sigma * i * n * a * x

    def total_cross_section(self, sigma: List[float]) -> float:
        """Returns the total cross-section, given the cross-sections of individual interactions"""
        return sum(sigma)

    def absorption_cross_section(self, sigma: List[float]) -> float:
        """Computes the absorption cross-section, given the cross-sections of individual
        interactions"""
        return sum(sigma)

    def total_cross_section_alt(self, sigma_s: float, sigma_a: float) -> float:
        """Alternate form of the total_cross_section formula. Calculates the total cross-section
        given the scattering cross-section and absorption cross-section"""
        return sigma_s + sigma_a

    def scattering_cross_section(self, sigma_e: float, sigma_i: float) -> float:
        """Returns the scattering cross-section given the elastic and inelastic scattering
        cross-sections"""
        return sigma_e + sigma_i

    def collision_density(self, i: float, n: float, sigma_t: float) -> float:
        """Computes the collision density given the beam intensity, target atom density, and total
        cross-section"""
        return i * n * sigma_t

    def collision_density_alt(self, i: float, sigma_macro: float) -> float:
        """An alternate version of the formula above, where the collision density is returned given
        the beam intensity and macroscopic cross-section"""
        return i * sigma_macro

    def macroscopic_cross_section(self, n: float, sigma: float) -> float:
        """Calculates the macroscopic cross-section of any neutron interaction given the atom
        density of the target, and the cross-section of the reaction/interaction type"""
        return n * sigma

    def intensity_thick_target(self, i_0: float, sigma_macro: float, x: float) -> float:
        """Returns the beam intensity at distance x into a thick target given the initial intensity
        and macroscopic cross-section. This can also be used to compute the intensity at the
        detector, where x is the distance to the detector from the beam source"""
        return i_0 * math.exp(-sigma_macro * x)

    def mean_free_path(self, sigma_macro: float) -> float:
        """Computes the neutron mean free path given the total macroscopic cross-section"""
        return 1 / sigma_macro

    def macroscopic_cross_section_mixture(
        self, n_x: float, sigma_x: float, n_y: float, sigma_y: float
    ) -> float:
        """Calculates the macroscopic cross-section of a homogeneous mixture of two nuclear
        materials given their atom densities and the cross-sections of both materials based on the
        nature of their interactions"""
        return n_x * sigma_x + n_y * sigma_y

    def macroscopic_cross_section_mixture_alt(
        self, sigma_macro_x: float, sigma_macro_y: float
    ) -> float:
        """Alternate version of the formula above, where the input args are the macroscopic
        cross-sections of each element"""
        return sigma_macro_x + sigma_macro_y

    def molecule_cross_section(
        self, m: int, sigma_m: float, n: int, sigma_n: float
    ) -> float:
        """Returns the cross-section of molecule X_m Y_n given the cross-sections of the particular
        interaction between each element"""
        return m * sigma_m + n * sigma_n

    def molecule_cross_section_alt(self, sigma_macro: float, n: float) -> float:
        """Alternate form of the equation above, given the macroscopic cross-section of the molecule
        and the number of molecules"""
        return sigma_macro / n

    def mb_collision_rate(self, macro_sigma_t: float, i: List[float]) -> float:
        """Returns the total interaction rate in a multi-beam apparatus given the total macroscopic
        cross-section and the intensities of all beams.
        """
        return macro_sigma_t * sum(i)

    def mb_collision_rate_densities(
        self, macro_sigma_t: float, v: float, n: Union[float, List[float]]
    ) -> float:
        """Calculates the total interaction rate in a multi-beam apparatus given the total
        macroscopic cross-section, the speed of the neutrons, and the densities of each neutron
        beams. In the case of an average or single density value, n can be a single input
        """
        if isinstance(list, n):
            return macro_sigma_t * v * sum(n)
        return macro_sigma_t * v * n

    def neutron_flux(self, n: float, v: float) -> float:
        """Computes the neutron flux given the density of the beam and the speed of the neutrons"""
        return n * v

    def collision_density_flux(self, sigma_macro_t: float, phi: float) -> float:
        """Returns the collision density given the total macroscopic cross-section and the neutron
        flux"""
        return sigma_macro_t * phi

    def potential_elastic_scattering(self, r: float) -> float:
        """Calculates the scattering cross-section due to forces exerted on the passing neutron from
        the target nucleus, or potential elastic scattering given the radius of the target nucleus
        """
        return 4 * math.pi * r**3

    def breit_wigner_one_level(
        self,
        gamma_r: float,
        g: float,
        gamma_n: float,
        gamma_g: float,
        e: float,
        e_r: float,
    ) -> float:
        """Computes the radiative capture cross-section in above the 1/v region, given the neutron
        wavelength, statistical factor, neutron width, radiation width, current energy level, and
        the energy of the neutrons"""
        return (gamma_r * g / (4 * math.pi)) * (
            gamma_n * gamma_g / ((e - e_r) ** 2 + ((gamma_n + gamma_g) ** 2) / 4)
        )

    def total_cross_section_resonances(self, r: float, c: float, e: float) -> float:
        """Calculates the total cross-section given the nuclear radius, constant C, and the neutron
        energy. The first term arises due to elastic scattering"""
        return self.potential_elastic_scattering(r) + c / math.sqrt(e)

    def energy_after_collision(self, e: float, a: int, theta: float) -> float:
        """Returns the post-collision energy given the initial energy, atomic mass number, and
        scattering angle, in radians"""
        return (
            (math.cos(theta) + math.sqrt(a**2 - math.sin(theta) ** 2)) ** 2
            * e
            / ((a + 1) ** 2)
        )

    def collision_parameter(self, a: int) -> float:
        """Computes the collision parameter given the atomic mass number"""
        return ((a - 1) / (a + 1)) ** 2

    def minimum_energy_after_collision(self, a: int, e: float) -> float:
        """Calculates the minimum possible energy after a collision, given the atomic mass number
        and energy before the collision"""
        return self.collision_parameter(a) * e

    def average_energy_scattering(self, alpha: float, e: float) -> float:
        """Returns the average energy of a scattered neutron given the collision parameter and
        initial neutron energy"""
        return (1 + alpha) * e / 2

    def average_energy_loss_scattering(self, alpha: float, e: float) -> float:
        """Computes the average energy lost by a neutron after scattering given the collision
        parameter, and initial energy"""
        return (1 - alpha) * e / 2

    def average_fractional_energy_loss(self, alpha: float) -> float:
        """Calculates the average fractional energy loss of a scattered neutron, given the collision
        parameter"""
        return (1 - alpha) / 2

    def lethargy(self, e_m: float, e: float) -> float:
        """Returns a neutron's lethargy given the energy of the highest neutron in the system and
        the energy of the neutron"""
        return math.log(e_m / e)

    def average_change_in_lethargy(self, a: int) -> float:
        """Computes the average change in lethargy of a neutron given the mass number of the target
        nucleus"""
        return 1 - ((a - 1) ** 2) * math.log((a + 1) / (a - 1)) / (2 * a)

    def absorption_rate(self, sigma_a_e_0: float, v: float, n: float) -> float:
        """Calculates the reaction rate for an absorption reaction given the macroscopic absorption
        cross-section at 0.0253 eV, neutron speed, and thermal neutron density rate"""
        return sigma_a_e_0 * self.neutron_flux(n, v)

    def absorption_rate_alt(self, sigma_a_e_0: float, phi: float) -> float:
        """Returns the reaction rate for an absorption reaction given the macroscopic absorption
        cross-section at 0.0253 eV, and the neutron flux"""
        return sigma_a_e_0 * phi

    def absorption_rate_non_1_v(
        self, g: float, sigma_a_e_0: float, v: float, n: float
    ) -> float:
        """Computes the absorption rate of a non-1/v absorber given the non-1/v factor, macroscopic
        absorption cross-section at 0.0253 eV, the speed of the neutrons, and the neutron density
        """
        return g * self.absorption_rate(sigma_a_e_0, v, n)

    def absorption_rate_non_1_v_alt(
        self, g: float, sigma_a_e_0: float, phi: float
    ) -> float:
        """Calculates the absorption rate of a non-1/v absorber given the non-1/v factor,
        macroscopic absorption cross-section at 0.0253 eV, and the neutron flux"""
        return g * self.absorption_rate_alt(sigma_a_e_0, phi)

    def capture_to_fission_ratio(self, sigma_gamma: float, sigma_f: float) -> float:
        """Returns the capture-to-fission ratio given the radiative capture and fission
        cross-sections"""
        return sigma_gamma / sigma_f

    def beta_ray_emission_rate(self, t: float) -> float:
        """Computes the rate of beta-ray emission from an interval of 10s to several weeks after
        fission occurs, given the time after fission in days. Returns the rate in beta-rays/s
        """
        return 3.8e-6 * t ** (-1.2)

    def gamma_ray_emission_rate(self, t: float) -> float:
        """Calculates the rate of gamma-ray emission from an interval of 10s to several weeks after
        fission occurs, given the time after fission in days. Retuns the rate in gamma-rays/s
        """
        return 1.9e-6 * t ** (-1.2)

    def fission_product_activity(self, t: float) -> float:
        """Returns the activity of fission products in Ci given the time in days after fission
        occurs"""
        return self.beta_ray_emission_rate(t) / conversions.activity["Bq/Ci"]

    def fission_product_activity_reactor(
        self, p: float, t_1: float, t_2: float
    ) -> float:
        """Computes the fission product activity in Ci for the following two scenarios:
        - The entire reactor: p is the (constant) operating power in MW, t_1 is the number of days
        that the reactor operated, and t_2 is the number of days since the reactor shutdown. The
        returned quantity is the fission product activity t_2 days after the reactor shutdown.
        - A single fuel rod: p is the power produced by the rod in MW, t_1 is the number days that
        the rod was in the reactor before removal, and t_2 is the days since removal. The returned
        quantity is the fission product activity of the rod, t_2 days after removal."""
        return 1.4e6 * p * (t_2 ** (-0.2) - (t_1 + t_2) ** (-0.2))

    def decay_energy_rate(self, t: float) -> float:
        """Calculates the decay energy rate of fission products in MeV/s given the time in days
        after the fission occured"""
        return 2.8e-6 * t ** (-1.2)

    def neutrons_released_fission_v1(self, nu: float, sigma_f: float, sigma_a) -> float:
        """Returns the number of neutrons released in fission per neutron absorbed by a fissile
        nucleus, given the average number of neutrons released per fission, fission cross-section,
        and absorption cross-section"""
        return nu * sigma_f / sigma_a

    def neutrons_released_fission_v2(
        self, nu: float, sigma_f: float, sigma_gamma: float
    ) -> float:
        """Alternate version of the above formula. Computes the number of neutrons released in
        fission per neutron absorbed by a fissile nucleus, given the average number of neutrons
        released per fission, fission cross-section, and radiative capture cross-section
        """
        return nu * sigma_f / (sigma_f + sigma_gamma)

    def neutrons_released_fission_v3(self, nu: float, alpha: float) -> float:
        """Alternate version of the above formula. Calculates the number of neutrons released in
        fission per neutron absorbed by a fissile nucleus, given the average number of neutrons
        released per fission and the capture-to-fission ratio"""
        return nu / (1 + alpha)

    def neutrons_released_fission_mixture(
        self, macro_sigma_a: float, nu: List[float], sigma_macro_f: List[float]
    ) -> float:
        """Returns the average number of neutrons emitted per neutron absorbed in a mixture of
        nuclides, given mixture's macroscopic absorption cross-section, average number of neutrons
        released per fission and macroscopic fission cross-sections of each nuclide"""
        return sum(nu[i] * sigma_macro_f[i] for i in range(len(nu))) / macro_sigma_a

    def prompt_neutron_spectrum(self, e: float) -> float:
        """Computes the probability or the dependence of the fraction of neutrons per MeV on neutron
        energy, given the neutron energy in MeV"""
        return 0.453 * math.exp(-1.036 * e) * math.sinh(2.29 * e)

    def fission_rate_235u_reactor(self, p: float) -> float:
        """Calculates the fissions per day an entire reactor operating at a thermal power of p MW
        containing uranium-235"""
        return 2.7e21 * p

    def burnup_rate_235u_reactor(self, p: float) -> float:
        """Returns the number of grams of uranium-235 fissioned per day in an entire reactor
        operating at a thermal power of p MW"""
        return (
            self.fission_rate_235u_reactor(p)
            * data.atomic_masses["235U"]
            / data.constants["N_A (1/(g-mol))"]
        )

    def scattered_photon_energy(self, e: float, theta: float) -> float:
        """Computes the energy of a photon scattered by the Compton effect in MeV, given the
        incident photon energy in MeV and the scattering angle in radians"""
        e_e = (
            data.constants["m_e (kg)"] * data.constants["c (m/s)"] ** 2
        ) / conversions.energy["J/MeV"]
        return e * e_e / (e * (1 - math.cos(theta)) + e_e)

    def scattered_photon_wavelength(self, lambda_: float, theta: float) -> float:
        """Calculates the wavelength of a photon scattered due to the Compton effect in cm, given
        the incident photon wavelength and angle of the scattered photon in radians"""
        return data.constants["lambda_C (cm)"] * (1 - math.cos(theta)) + lambda_

    def compton_cross_section_per_atom(self, z: int, e_sigma_c: float) -> float:
        """Returns the Compton cross-section per atom, given the number of electrons per atom and
        the Compton cross-section per electron"""
        return z * e_sigma_c

    def gamma_ray_cross_section(
        self, sigma_pe: float, sigma_pp: float, sigma_c: float
    ) -> float:
        """Computes the total gamma-ray interaction cross-section given the photoelectric effect,
        pair production, and Compton effect cross-sections"""
        return sigma_pe + sigma_pp + sigma_c

    def attenuation_coefficient_v1(self, n: float, sigma: float) -> float:
        """Calculates the attenuation coefficient given the atom density and the total cross-section
        of gamma-ray interactions. Calculated similarly to neutron cross-sections."""
        return self.macroscopic_cross_section(n, sigma)

    def attenuation_coefficient_v2(
        self, n: float, sigma_pe: float, sigma_pp: float, sigma_c: float
    ) -> float:
        """Returns the attenuation coefficient given the atom density, the photoelectric effect,
        pair production, and Compton effect cross-sections"""
        sigma = self.gamma_ray_cross_section(sigma_pe, sigma_pp, sigma_c)
        return self.macroscopic_cross_section(n, sigma)

    def mass_attenuation_coefficient(
        self, mu_pe: float, mu_pp: float, mu_c: float, rho: float
    ) -> float:
        """Computes the mass attenuation coefficient given the the photoelectric effect, pair
        production, and Compton effect attenuation coefficients, and the physical density
        """
        return (mu_pe + mu_pp + mu_c) / rho

    def mass_attenuation_coeff_compton(
        self, a: int, m: float, e_sigma_c: float
    ) -> float:
        """Calculates the mass attenuation coefficient in energy ranges where Compton scattering
        is the dominant mode of interaction, given the atomic mass number, gram atomic weight, and
        Compton scattering cross-section per electron"""
        return data.constants["N_A (1/(g-mol))"] * a * e_sigma_c / m

    def attenuation_coefficient_mixture(self, mu: List[float]) -> float:
        """Returns the attenuation coefficient of a mixture of elements, given the attenuation
        coefficents of each element"""
        return sum(mu)

    def mass_attenuation_coeff_mixture(
        self, omega: List[float], mu_star: List[float]
    ) -> float:
        """Computes the mass attenuation coefficient of a mixture given the weight percents and mass
        attenuation coefficients of each element"""
        return sum(mu_star[i] * omega[i] for i in range(len(mu_star))) / 100

    def gamma_ray_mean_free_path(self, mu: float) -> float:
        """Calculates the mean free path of a gamma ray given the attenuation coefficient"""
        return 1 / mu

    def gamma_ray_intensity_v1(self, i_0: float, mu: float, x: float) -> float:
        """Returns the intensity of the photons that penetrate a target without having a collision,
        given the mono-energetic gamma-ray beam intensity, attenuation coefficient, and the
        thickness of the target"""
        return i_0 * math.exp(-mu * x)

    def gamma_ray_intensity_v2(
        self, i_0: float, mu_star: float, rho: float, x: float
    ) -> float:
        """Alternate form of the equation above. Computes the intensity of the photons that
        penetrate a target without having a collision, given the mono-energetic gamma-ray beam
        intensity, mass attenuation coefficient, the physical density of the target, and the
        thickness of the target"""
        return i_0 * math.exp(mu_star * rho * x)

    def gamma_ray_collision_density(self, i: float, mu: float) -> float:
        """Calculates the collision density of a gamma-ray given its intensity and total attenuation
        coefficient"""
        return i * mu

    def avg_energy_deposited_compton(
        self, t_bar: float, i: float, mu_c: float
    ) -> float:
        """Returns the average energy deposited by Compton scattering, given the average energy of
        the recoiling electron, gamma-ray intensity, and Compton attenuation coefficient
        """
        return t_bar * i * mu_c

    def compton_absorption_cross_section(
        self, e: float, t_bar: float, sigma_c: float
    ) -> float:
        """Computes the Compton absorption cross-section given the energy of the gamma-rays,
        average energy of the recoiling electron, and the Compton scattering cross-section
        """
        return t_bar * sigma_c / e

    def compton_absorption_coefficient(
        self, e: float, t_bar: float, mu_c: float
    ) -> float:
        """Calculates the Compton absorption coefficient given the energy of the gamma-rays, average
        energy of the recoiling electron, and the Compton attenuation coefficient"""
        return t_bar * mu_c / e

    def compton_energy_deposition_rate(self, e: float, i: float, mu_ca: float) -> float:
        """Returns the energy deposition rate per unit volume by Compton scattering given the
        gamma-ray energy, intensity, and the Compton absorption coefficient"""
        return e * i * mu_ca

    def linear_absorption_coefficient(
        self, mu_pe: float, mu_pp: float, mu_ca: float
    ) -> float:
        """Computes the linear absorption coefficient given the attenuation coefficients of the
        photelectric effect, pair production, and the Compton absorption coefficient"""
        return mu_pe + mu_pp + mu_ca

    def total_energy_deposition_rate_v1(
        self, e: float, i: float, mu_pe: float, mu_pp: float, mu_ca
    ) -> float:
        """Calculates the total energy deposition rate per unit volume given the gamma-ray energy,
        intensity, and linear absorption coefficient"""
        return e * i * self.linear_absorption_coefficient(mu_pe, mu_pp, mu_ca)

    def total_energy_deposition_rate_v2(self, e: float, i: float, mu_a: float) -> float:
        """Returns the total energy deposition rate per unit volume given the gamma-ray energy,
        intensity, the attenuation coefficients of the photelectric effect, pair production, and the
        linear absorption coefficient"""
        return e * i * mu_a

    def energy_deposition_rate_mass(
        self, e: float, i: float, mu_a: float, rho: float
    ) -> float:
        """Computes the energy deposition rate per unit mass given the energy and intensity of the
        gamma-rays, the linear absorption coefficient, and the physical density"""
        return e * i * mu_a / rho

    def energy_deposition_rate_att(self, e: float, i: float, mu_rho_a: float) -> float:
        """Calculates the energy deposition rate per unit mass given the energy and intensity of the
        gamma-rays, and the mass absorption coefficient"""
        return e * i * mu_rho_a

    def stopping_power(self, de_dx_col: float, de_dx_rad: float) -> float:
        """Returns the stopping power if nuclear reactions between charged particles do not occur,
        given the energy loss per unit distance due to collisions and the energy loss by radiation
        """
        return de_dx_col + de_dx_rad

    def bragg_kleeman_rule(self, m: float, rho: float, r_a: float) -> float:
        """Computes the range of an alpha-particle and other materials in cm given the atomic mass,
        physical density, and the range of particles in air"""
        return 4.2e-4 * math.sqrt(m) * r_a / rho

    def brag_kleeman_rule_mixture(
        self, m: List[float], gamma: List[float], rho: float, r_a: float
    ) -> float:
        """Alternate form of the Bragg-Kleeman rule for mixtures, where the range of a mixture of
        elements is given in cm. Requires a list of the atomic masses, atom fractions (summing to
        1), the average density, and the range of the air medium"""
        m = sum(m[i] * gamma[i] for i in range(len(m))) ** 2
        return self.bragg_kleeman_rule(m, rho, r_a)

    def relative_stopping_power_v1(self, r: float, r_a: float) -> float:
        """Calculates the relative stopping power (unitless) given the range of the particles, and
        the range of particles in air"""
        return r_a / r

    def relative_stopping_power_v2(self, rho: float, m: float) -> float:
        """Alternate form of the relative stopping power expression, given the physical density and
        atomic mass of the particle"""
        return 3100 * rho / math.sqrt(m)

    def specific_ionization_v1(self, i_0: float, mu: float, x: float) -> float:
        """Returns the specific ionization of a beam of beta-rays given the initial ionization
        (at x = 0), the attenuation coefficient, and distance into the absorber"""
        return i_0 * math.exp(-mu * x)

    def specific_ionization_v2(
        self, i_0: float, mu_star: float, rho: float, x: float
    ) -> float:
        """Alternate form of the specific ionization equation. Computes the specific ionization of
        a beta-ray beam given the initial ionization, absorption coefficient, density of the
        medium, and the distance into the absorber"""
        return i_0 * math.exp(-mu_star * rho * x)

    def beta_ray_attenuation_coeff(self, e_max: float) -> float:
        """Calculates the mass attenuation coefficient of beta-rays given the max beta-ray energy,
        in MeV"""
        return 17 / (e_max**1.14)

    def beta_ray_max_range(self, e_max: float, rho: float) -> float:
        """An empirical formula that returns the maximum range in cm of beta-rays given the maximum
        energy of the beta-rays, and the density of the medium"""
        return (
            0.412 * (e_max ** (1.265 - 0.0954 * math.log(e_max))) / rho
            if e_max < 2.5
            else (0.530 * e_max - 0.106) / rho
        )
