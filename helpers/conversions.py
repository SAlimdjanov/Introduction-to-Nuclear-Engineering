"""
conversions.py

Notes:
- All factors are from Appendix 1 unless stated otherwise
- The L/mol conversion factor is only applicable to an ideal gas
- Bq/Ci is equivalent to (disintegrations/s)/Ci
- The power_density and energy_flux tables can also be taken as heat source density and power flux,
  respectively.

"""

from math import pi


length = {
    "cm/m": 100,
    "km/m": 1e-3,
    "in/m": 39.37,
    "ft/m": 3.281,
    "mi/m": 6.214e-4,
}

area = {
    "cm^2/m^2": 1e4,
    "km^2/m^2": 1e-6,
    "in^2/m^2": 1550,
    "ft^2/m^2": 10.76,
    "mi^2/m^2": 3.861e-7,
    "b/cm^2": 1e24,
}

volume = {
    "cm^3/m^3": 1e6,
    "L/m^3": 1e3,
    "in^3/m^3": 6.102e4,
    "ft^3/m^3": 35.31,
    "L/mol": 22.4,
}

mass = {
    "g/kg": 1e3,
    "lb/kg": 2.205,
    "t(short)/kg": 1.102e-3,
    "t(metric)/kg": 1e-3,
    "g/amu": 1.6606e-24,
    "g/electron": 9.1095e-28,
}

time = {
    "s/min": 60,
    "s/hr": 3600,
    "s/day": 24 * 3600,
    "s/mo": 2.629746e6,
    "s/yr": 365 * 24 * 3600,
}

energy = {
    "erg/J": 1e7,
    "(g-cal)/J": 0.2388,
    "Btu/J": 9.418e-4,
    "J/kWh": 3.6e6,
    "J/eV": 1.6022e-19,
    "J/MeV": 1.6022e-13,
    "MeV/amu": 931.502,
}

power = {
    "kW/W": 1e-3,
    "MW/W": 1e-6,
    "(Btu/hr)/W": 3.412,
    "(MeV/s)/W": 6.242e12,
}

power_density = {
    "(W/cm^3)/(W/m^3)": 1e-6,
    "[cal/(cm^3-s)]/(W/m^3)": 2.388e-7,
    "[Btu/(ft^3-hr)]/(W/m^3)": 0.09662,
    "[MeV/(cm^3-s)]/(W/m^3)": 6.242e6,
}

energy_flux = {
    "(W/cm^2)/(W/m^2)": 1e-4,
    "[cal/(cm^2-s)]/(W/m^2)": 2.388e-5,
    "[Btu/(ft^2-hr)]/(W/m^2)": 0.3170,
    "[MeV/(cm^2-s)]/(W/m^2)": 6.242e8,
}

thermal_conductivity = {
    "[W/(cm-C)]/[W/(m-K)]": 0.01,
    "[Btu/(ft-hr-F)]/[W/(m-K)]": 2.388e-3,
    "[(Btu-ft)/(hr-ft^2-F)]/[W/(m-K)]": 0.5778,
}

viscosity = {
    "[g/(cm-s)]/[kg/(m-s)]": 10,  # This quantity is also known as Poise, denoted by P
    "cP/[kg/(m-s)]": 1e3,  # Centipoise
    "[lb/(ft-hr)]/[kg/(m-s)]": 2419,
}

pressure = {
    "(dyne/cm^2)/Pa": 10,
    "psi/Pa": 1.45e-4,
    "bar/Pa": 1e-5,
    "Pa/atm": 101325,
}

activity = {
    "Bq/Ci": 3.7e10,
}

unit_prefixes = {
    "P": 1e15,
    "T": 1e12,
    "G": 1e9,
    "M": 1e6,
    "k": 1e3,
    "c": 1e-2,
    "m": 1e-3,
    "mu": 1e-6,
    "n": 1e-9,
    "p": 1e-12,
    "f": 1e-15,
}

angle = {
    "deg/rad": 180 / pi,
}


def temperature(t: float, input_units="K", output_units="C") -> float:
    """Converts an inpit temperature to degrees Celsius, Kelvin, or Fahrenheit"""
    units = input_units, output_units
    # Convert the input temperature to Celsius as an intermediate step
    match units:
        case "C", "K":
            temp = t + 273.15
        case "K", "C":
            temp = t - 273.15
        case "C", "F":
            temp = (t * 9 / 5) + 32
        case "F", "C":
            temp = (t - 32) * 5 / 9
        case "F", "K":
            temp = (t - 32) * 5 / 9 + 273.15
        case "K", "F":
            temp = (t - 273.15) * 9 / 5 + 32
        case _:
            raise ValueError(
                "Units must be 'C' (Celsius), 'K' (Kelvin), or 'F' (Fahrenheit)"
            )
    return temp
