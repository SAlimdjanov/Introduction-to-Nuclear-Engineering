"""
conversions.py

Notes:
- The L/mol conversion factor is only applicable to an ideal gas
- Bq/Ci is equivalent to (disintegrations/s)/Ci

"""

energy = {
    "J/eV": 1.6022e-19,
    "J/MeV": 1.6022e-13,
    "eV/MeV": 1e6,
    "Btu/J": 9.418e-4,
    "J/kWh": 3.6e6,
    "MeV/amu": 931.502,
}

mass = {
    "g/amu": 1.6606e-24,
    "g/electron": 9.1095e-28,
    "g/kg": 1000,
}

volume = {
    "L/mol": 22.4,
    "mL/L": 1e3,
    "cm^3/m^3": 1e6,
}

time = {
    "s/us": 1e-6,
    "s/min": 60,
    "s/hr": 3600,
    "s/day": 24 * 3600,
    "s/mo": 2.629746e6,
    "s/yr": 365 * 24 * 3600,
}

length = {
    "m/fm": 1e-15,
    "m/km": 1000,
    "m/um": 1e-6,
}

pressure = {
    "Pa/atm": 101325,
    "bar/Pa": 1e-5,
}

activity = {
    "Bq/Ci": 3.7e10,
}


def temperature(t: float, units="C") -> float:
    """Convert temperature to and from degrees Celsius and degrees kelvin. Kwarg corresponds to the
    input unit of temperature"""
    if units not in ("C", "K"):
        raise ValueError("Units must be 'C' (Celsius) or 'K' (Kelvin)")
    diff = 273.15
    return t + diff if units == "C" else t - diff
