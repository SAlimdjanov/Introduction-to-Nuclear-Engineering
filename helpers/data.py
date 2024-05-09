"""
data.py

"""

constants = {
    "amu (g)": 1.66054e-24,
    "amu (MeV)": 931.494,
    "N_A (1/(g-mole))": 0.6022137e24,
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
    "2H": 0.015,
    "16O": 99.759,
    "17O": 0.037,
    "18O": 0.204,
}

atomic_weights = {  # Unitless
    "16O": 15.99492,
    "17O": 16.99913,
    "18O": 17.99916,
    "Na": 22.98977,
    "Cl": 35.45270,
    "H": 1.00797,
    "O": 15.9994,
    "235U": 235.0439,
    "238U": 238.0508,
}

neutral_atomic_masses = {  # amu
    "n": 1.008665,
    "1H": 1.007825,
    "2H": 2.014102,
    "3H": 3.016049,
    "4He": 4.002604,
    "12C": 12.000000,
    "13C": 13.003354,
    "14C": 14.003242,
}

# Alpha particle spectrum of 226Ra. Keys are alpha particle energies and values are the relative
# number of particles, expressed as a percentage.
alpha_spectrum_226ra = {
    "4.782": 94.6,
    "4.599": 5.4,
    "4.340": 5.1e-3,
    "4.194": 7e-4,
}
