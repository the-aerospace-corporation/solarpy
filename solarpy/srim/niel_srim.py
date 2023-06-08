import numpy as np
import scipy as sp
from scipy import integrate


def get_energy_vacancy(T_d):
    """
    Calculates the energy necessary to produce a vacancy

    Args:
        T_d: threshold energy for displacement in eV of target material (10eV for GaAs)

    Returns:
        Energy necessary to produce a vacancy in keV/Vacancy
    """
    energy_vacancy = (1/1000) * ((T_d/0.4) +2)

    return energy_vacancy


def get_NIEL_depth(energy_vacancy, ion_v, rec_v, density):
    niel_depth = np.zeros((len(ion_v[:,0]),2))
    niel_depth[:,0] = ion_v[:,0]
    niel_depth[:,1] = energy_vacancy*(ion_v[:,1] + rec_v[:,1]) *(1e5/density)

    return niel_depth

def get_energy_lost(energy_vacancy, ion_v, rec_v, ion_i, rec_i):
    energy_lost = np.zeros((len(ion_v[:, 0]), 2))
    energy_lost[:, 0] = ion_v[:, 0]
    energy_lost[:, 1] = energy_vacancy * (ion_v[:,1] + rec_v[:,1]) + 1e-3*(ion_i[:,1] + rec_i[:,1])
    return energy_lost

def get_cumulative_energy_loss(energy_lost):
    Ec_D = sp.integrate.cumtrapz(energy_lost[:, 1], energy_lost[:, 0], initial=0)
    return Ec_D

def get_niel_from_srim(particle_energy_MeV, T_d, ion_v, rec_v, ion_i, rec_i, density):
    energy_vacancy = get_energy_vacancy(T_d)
    niel_depth = get_NIEL_depth(energy_vacancy, ion_v, rec_v, density)
    energy_lost = get_energy_lost(energy_vacancy, ion_v, rec_v, ion_i, rec_i)
    Ec_D = get_cumulative_energy_loss(energy_lost)
    niel_srim = np.zeros((len(Ec_D),2))
    niel_srim[:,0] = particle_energy_MeV-Ec_D/1000
    niel_srim[:,1] = niel_depth[:,1]
    return niel_srim
