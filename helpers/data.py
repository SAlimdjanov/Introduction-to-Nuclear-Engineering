"""
data.py

Notes:
- All values were obtained from the book unless otherwise stated
- Units associated with quantities are commented at the top of each dictionary, where appropriate
- Avagadro's number can also be expressed in units of atoms/mol
- The atomic_masses and half_lives dictionaries contain a combination of values from Appendix II
  and the Table of Nuclides found at https://www-nds.iaea.org/relnsd/vcharthtml/VChartHTML.html

"""

from .conversions import conversion_factors


constants = {
    "amu (g)": 1.66054e-24,
    "amu (MeV)": 931.494,
    "N_A (1/(g-mol))": 0.6022137e24,
    "k_B (J/K)": 1.38066e-23,
    "lambda_C (eV/K)": 8.61707e-5,
    "lambda_C (cm)": 2.42631e-10,
    "m_e (kg)": 9.10939e-31,
    "m_e (amu)": 5.485799e-4,
    "m_e (MeV)": 0.510999,
    "q_e (C)": 1.602192e-19,
    "m_n (kg)": 1.674929e-27,
    "m_n (amu)": 1.008665,
    "m_n (MeV)": 939.56563,
    "h (J-s)": 6.626075e-34,
    "h (eV-s)": 4.13572e-15,
    "m_p (kg)": 1.67262e-27,
    "m_p (amu)": 1.007276,
    "m_p (MeV)": 938.27231,
    "c (m/s)": 2.997925e8,
    "alpha (MeV)": 15.56,
    "beta (MeV)": 17.23,
    "gamma (MeV)": 0.697,
    "zeta (MeV)": 23.285,
    "delta (MeV)": 12.0,
}

abundances = {  # atom percent, a/o
    "1H": 99.986,
    "2H": 0.015,
    "12C": 100.0,
    "16O": 99.759,
    "17O": 0.037,
    "18O": 0.204,
    "234U": 0.0057,
    "235U": 0.72,
    "238U": 99.27,
}

atomic_masses = {  # amu
    "n": 1.00867,
    "p": 1.00728,
    "gamma": 0.00000,
    "alpha": 4.00151,
    "H": 1.00797,
    "C": 12.01070,
    "O": 15.99940,
    "Na": 22.98977,
    "Cl": 35.45270,
    "1H": 1.00783,
    "2H": 2.01410,
    "3H": 3.01605,
    "3He": 3.01603,
    "4He": 4.00260,
    "5Li": 5.01254,
    "6Li": 6.01512,
    "7Li": 7.01600,
    "8Be": 8.00535,
    "9Be": 9.01218,
    "12C": 12.00000,
    "13C": 13.00335,
    "14C": 14.00324,
    "14N": 14.00307,
    "16O": 15.99492,
    "17O": 16.99913,
    "18O": 17.99916,
    "50V": 49.94716,
    "51V": 50.94396,
    "59Co": 58.93319,
    "60Co": 59.93382,
    "90Sr": 89.90774,
    "115In": 114.90388,
    "116In": 115.90526,
    "206Pb": 205.97447,
    "207Pb": 206.97590,
    "208Pb": 207.97665,
    "138Ba": 137.90525,
    "222Ra": 222.01537,
    "226Ra": 226.02541,
    "234U": 234.04095,
    "235U": 235.04393,
    "238U": 238.05079,
    "239Pu": 239.05216,
    "240Pu": 240.05381,
    "241Pu": 241.05685,
    "242Pu": 242.05874,
}

alpha_spectrum_226ra = {  # "MeV": %
    "4.782": 94.6,
    "4.599": 5.4,
    "4.340": 5.1e-3,
    "4.194": 7e-4,
}

half_lives = {  # seconds
    "n": 12 * conversion_factors["s/min"],
    "3H": 12.33 * conversion_factors["s/yr"],
    "14C": 5736 * conversion_factors["s/yr"],
    "60Co": 1925.28 * conversion_factors["s/day"],
    "90Sr": 28.91 * conversion_factors["s/yr"],
    "134Cs": 2.0652 * conversion_factors["s/yr"],
    "135I": 6.7 * conversion_factors["s/hr"],
    "135Xe": 9.17 * conversion_factors["s/hr"],
    "137Cs": 30.08 * conversion_factors["s/yr"],
    "143Pm": 53.1 * conversion_factors["s/hr"],
    "214Bi": 19.71 * conversion_factors["s/min"],
    "214Pb": 27.06 * conversion_factors["s/min"],
    "214Po": 163.46 * conversion_factors["s/us"],
    "218Po": 3.097 * conversion_factors["s/min"],
    "222Rn": 3.8222 * conversion_factors["s/day"],
    "226Ra": 1600 * conversion_factors["s/yr"],
    "232Th": 1.41e10 * conversion_factors["s/yr"],
    "233Th": 23.3 * conversion_factors["s/min"],
    "233U": 1.592e5 * conversion_factors["s/yr"],
    "234U": 2.46e5 * conversion_factors["s/yr"],
    "235U": 7.038e8 * conversion_factors["s/yr"],
    "236U": 2.34e7 * conversion_factors["s/yr"],
    "238U": 4.68e9 * conversion_factors["s/yr"],
    "239U": 23.5 * conversion_factors["s/min"],
    "239Pu": 24110 * conversion_factors["s/yr"],
    "240Pu": 6564 * conversion_factors["s/yr"],
    "241Pu": 14.35 * conversion_factors["s/yr"],
    "242Pu": 3.733e5 * conversion_factors["s/yr"],
}
