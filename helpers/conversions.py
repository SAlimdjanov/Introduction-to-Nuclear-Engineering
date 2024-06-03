"""
conversions.py

Notes:
- The L/mol conversion factor is only applicable to an ideal gas
- Bq/Ci is equivalent to (disintegrations/s)/Ci
- Unit prefixes are included for easy conversions.

"""

energy = {
    "J/eV": 1.6022e-19,
    "J/MeV": 1.6022e-13,
    "Btu/J": 9.418e-4,
    "J/kWh": 3.6e6,
    "MeV/amu": 931.502,
}

mass = {
    "g/amu": 1.6606e-24,
    "g/electron": 9.1095e-28,
}

area = {
    "b/cm^2": 1e24,
}

volume = {
    "L/mol": 22.4,
    "cm^3/m^3": 1e6,
}

time = {
    "s/min": 60,
    "s/hr": 3600,
    "s/day": 24 * 3600,
    "s/mo": 2.629746e6,
    "s/yr": 365 * 24 * 3600,
}

pressure = {
    "Pa/atm": 101325,
    "bar/Pa": 1e-5,
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
