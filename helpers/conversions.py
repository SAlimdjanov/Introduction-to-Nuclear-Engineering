"""
conversions.py

"""

conversion_factors = {
    "J/MeV": 1.6022e-13,
    "g/amu": 1.6606e-24,
    "g/electron": 9.1095e-28,
    "MeV/amu": 931.502,
    "g/kg": 1000,
    "eV/MeV": 1e6,
    "Btu/J": 9.418e-4,
    "J/kWh": 3.6e6,
}


def convert_temperature(t: float, units="C") -> float:
    """Convert temperature to and from degrees Celsius and degrees kelvin. Kwarg corresponds to the
    input unit of temperature"""
    if units not in ("C", "K"):
        raise ValueError("Units must be 'C' (Celsius) or 'K' (Kelvin)")
    diff = 273.15
    return t + diff if units == "C" else t - diff
