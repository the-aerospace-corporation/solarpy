import numpy as np
import scipy as sp
from scipy import interpolate
from scipy import integrate


def get_relative_fluence_JPL_Aero(relative_damage_coefficients, differential_particle_spectrum):
    ###### DEPRECATED ######
    # Need to test if you need to use logFluence for differential electron spectrum and loglog of differential proton spectrum
    """
    Calculates the 1 MeV electron fluence using the JPL EQFLUX method using the differential particle spectrum. The original JPL EQFLUX program takes the integral particle spectrum and does what seems like a step-by-step differentiation then a Riemann sum to arrive at the total 1 MeV fluence.  This functions also does a Riemann sum but takes out the differentiation of the integral spectrum and simply uses the differential spectrum. The differential particle spectrum for an environment is easily obtained from Spenvis on AE8/AP8 and AE9/AP9 radiation environment models.  This function is similar to the 'oneMeVfluence_trapezoidalIntegration' function. The difference being this function uses Riemann sum integral appoximation and the 'oneMeVfluence_trapezoidalIntegration' function uses trapezoidal integration.

    Args:
        relative_damage_coefficients (ndarray) : 2d array where column 0 contains the particle energies and column 1 contains the relative damage coefficients
        differential_particle_spectrum (ndarray) : 2d array where column 0 contains the particle energies and column 1 contains the differential particle fluence

    Returns:
        Total 1MeV fluence

    """
    interpolate_relative_damage_coefficients = sp.interpolate.interp1d(np.log(relative_damage_coefficients[:, 0]),
                                                                       np.log(relative_damage_coefficients[:, 1]),
                                                                       fill_value='extrapolate')
    particle_energies = differential_particle_spectrum[:, 0]
    fluence = differential_particle_spectrum[:, 1]
    rdc = np.nan_to_num(np.exp(interpolate_relative_damage_coefficients(np.log(particle_energies))))
    oneMeVFluence = (fluence[:-1] * rdc[:-1] * ((np.exp(np.log(particle_energies[1:])) - np.exp(np.log(particle_energies[:-1]))))).sum()
    return oneMeVFluence


def get_relative_MeV_fluence_trapezoidal_integration(relative_damage_coefficients, differential_particle_spectrum):
    """
    Calculates the total 1MeV fluence given the relative damage coefficients and differential particle spectrum. The relative damage coefficients are interpolated for particle energies from the differential particle spectrum. The interpolation of the relative damage coefficients is accomplished by linear interpolation of the loglog of the relative damage coefficient. Relative damage coefficients outside the range of user supplied values is linearly extrapolated.  The interpolation and extrapolation of the relative damage coefficients is in accordance with the Solar Cell and GaAs Solar Cell Radiation Handbooks.  The interpolated relative damage coefficients are multiplied by the differential fluence to yield the 1MeV fluence at each particle energy.  The energy vs 1MeV fluence data is then integrated over all particle energies to yield the total 1MeV fluence

    Args:
        relative_damage_coefficients (ndarray) : 2d array where column 0 contains the particle energies and column 1 contains the relative damage coefficients
        differential_particle_spectrum (ndarray) : 2d array where column 0 contains the particle energies and column 1 contains the differential particle fluence

    Returns:
        Total 1MeV fluence

    """
    energy_vs_1MeV_fluence = get_energy_vs_relative_fluence(relative_damage_coefficients, differential_particle_spectrum)
    # oneMeV_fluence = sp.integrate.trapz(energy_vs_1MeV_fluence[:, 1], energy_vs_1MeV_fluence[:, 0])
    oneMeV_fluence = np.trapz(energy_vs_1MeV_fluence[:, 1], energy_vs_1MeV_fluence[:, 0])
    return oneMeV_fluence


def get_cumulative_relative_fluence(relative_damage_coefficients, differential_particle_spectrum, norm = True):
    """
    Calculates the cumulative 1MeV fluence for each particle energy

    Args:
        relative_damage_coefficients(2d ndarray): 2d array where column 0 contains the particle energies and column 1 contains the relative damage coefficients
        differential_particle_spectrum(2d ndarray): 2d array where column 0 contains the particle energies and column 2 contains the differential particle fluence

    Returns:
        2d nd array where column 0 is the particle energies and column 1 is the cumulative 1 Mev fluence.

    """
    energy_vs_1MeV_fluence = get_energy_vs_relative_fluence(relative_damage_coefficients, differential_particle_spectrum)
    cumulative_1MeV_fluence = np.zeros(np.shape(energy_vs_1MeV_fluence))[:-1]
    cumulative_1MeV_fluence[:, 0] = energy_vs_1MeV_fluence[:-1, 0]
    cumulative_1MeV_fluence[:, 1] = sp.integrate.cumtrapz(energy_vs_1MeV_fluence[:, 1], energy_vs_1MeV_fluence[:, 0])
    if norm:
        cumulative_1MeV_fluence[:, 1] = cumulative_1MeV_fluence[:, 1] / np.max(cumulative_1MeV_fluence[:, 1])
    return cumulative_1MeV_fluence


def get_energy_vs_relative_fluence(relative_damage_coefficients, differential_particle_spectrum):
    """
    Calculates the energy vs 1 MeV fluence data.  The relative damage coefficients are interpolated for particle energies from the differential particle spectrum. The interpolation of the relative damage coefficients is accomplished by linear interpolation of the loglog of the relative damage coefficient. Relative damage coefficients outside the range of user supplied values is linearly extrapolated.  The interpolation and extrapolation of the relative damage coefficients is in accordance with the Solar Cell and GaAs Solar Cell Radiation Handbooks. The interpolated relative damage coefficients are multiplied by the differential fluence to yield the 1MeV fluence at each particle energy.

    Args:
        relative_damage_coefficients(2d ndarray): 2d array where column 0 contains the particle energies and column 1 contains the relative damage coefficients
        differential_particle_spectrum(2d ndarray): 2d array where column 0 contains the particle energies and column 2 contains the differential particle fluence

    Returns:
        2d ndarray where column 0 is the particle energies and column 1 is the 1 MeV electron fluence. The particle energies are from the differential particle spectrum as the relative damage coefficients are interpolated and extropolated the differential particle spectrum.

    """
    # intp_relative_damage_coefficients = sp.interpolate.interp1d(np.log(relative_damage_coefficients[:, 0]), np.log(relative_damage_coefficients[:, 1]), fill_value='extrapolate')
    particle_energies = differential_particle_spectrum[:, 0]
    fluence = differential_particle_spectrum[:, 1]
    # rdc = np.exp(intp_relative_damage_coefficients(np.log(particle_energies)))
      # Convert those pesky NANs to 0s
    rdc = np.exp(np.interp(np.log(particle_energies), np.log(relative_damage_coefficients[:, 0]), np.log(relative_damage_coefficients[:, 1])))
    if np.isnan(rdc).any():
        rdc = np.nan_to_num(rdc)
    energy_vs_1MeV_fluence = np.vstack((particle_energies, rdc * fluence)).T
    return energy_vs_1MeV_fluence


def interpolate_electron_spectrum(integral_differential_electron_spectrum, energy_increment=0.01):
    """
    Interpolates the electron particle spectrum according to The GaAs Radiation Handbook.  The Handbook suggests interpolating the integral electron spectrum using linear interpolation of the "particle energy vs log(electron fluence)"

    Args:
        integral_differential_electron_spectrum (ndarray): 2d numpy array where column 0 is the particle energy and column 1 is the electron fluence
        energy_increment (float): spacing in which to linearly interpolate spectrum

    Returns:
        A 2d numpy array of the interpolated electron spectrum where column 0 is the particle energy and column 1 is the electron fluence

    """
    intperolate_electron_spectrum = sp.interpolate.interp1d(integral_differential_electron_spectrum[:, 0], np.log(integral_differential_electron_spectrum[:, 1]))
    particle_energy = np.arange(np.min(integral_differential_electron_spectrum[:, 0]), np.max(integral_differential_electron_spectrum[:, 0]), energy_increment)
    fluence = np.exp(intperolate_electron_spectrum(particle_energy))
    fluence = np.nan_to_num(fluence)
    intepolated_electron_spectrum = np.vstack((particle_energy, fluence)).T
    return intepolated_electron_spectrum


def interpolate_proton_spectrum(integral_differential_proton_spectrum, energy_increment=0.001):
    """
    Interpolates the proton particle spectrum according to The GaAs Radiation Handbook.  The Handbook suggests interpolating the integral pront spectrum using linear interpolation of the "log(particle energy) vs log(electron fluence)"

    Args:
        integral_differential_proton_spectrum (ndarray): 2d numpy array where column 0 is the particle energy and column 1 is the proton fluence
        energy_increment (float): spacing in which to linearly interpolate spectrum

    Returns:
        A 2d numpy array of the interpolated electron spectrum where column 0 is the particle energy and column 1 is the proton fluence

    """
    interpolate_proton_spectrum = sp.interpolate.interp1d(np.log(integral_differential_proton_spectrum[:, 0]), np.log(integral_differential_proton_spectrum[:, 1]))
    particle_energy = np.arange(np.min(np.log(integral_differential_proton_spectrum[:, 0])), np.max(np.log(integral_differential_proton_spectrum[:, 0])), energy_increment)
    fluence = np.exp(interpolate_proton_spectrum(particle_energy))
    fluence = np.nan_to_num(fluence)
    interpolated_proton_spectrum = np.vstack((np.exp(particle_energy), fluence)).T
    return interpolated_proton_spectrum
