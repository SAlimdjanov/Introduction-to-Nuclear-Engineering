"""
conversions.py

Notes:
- The L/mol conversion factor is only applicable to an ideal gas
- Bq/Ci is equivalent to (disintegrations/s)/Ci

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
    "Bq/Ci": 3.7e10,
    "L/mol": 22.4,
    "mL/L": 1e3,
    "Pa/atm": 101325,
    "s/us": 1e-6,
    "s/min": 60,
    "s/hr": 3600,
    "s/day": 24 * 3600,
    "s/mo": 2.629746e6,
    "s/yr": 365 * 24 * 3600,
}


def convert_temperature(t: float, units="C") -> float:
    """Convert temperature to and from degrees Celsius and degrees kelvin. Kwarg corresponds to the
    input unit of temperature"""
    if units not in ("C", "K"):
        raise ValueError("Units must be 'C' (Celsius) or 'K' (Kelvin)")
    diff = 273.15
    return t + diff if units == "C" else t - diff
