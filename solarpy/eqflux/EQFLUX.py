import numpy as np
import scipy as sp
from scipy import interpolate


def get_log_electron_integral_spectrum(electron_integral_particle_spectrum):
    """
    Generates the electron particle energy vs ln(fluence) spectrum of the integral particle spectrum in accordance with the JPL EQGAFLUX FORTRAN program as described in the GaAs Solar Cell Radiation Handbook

    Args:
        electron_integral_particle_spectrum (ndarray): 2d array of the integral electron particle spectrum where column 0 is the electron particle energy and column 1 is the integral fluence.  The integral particle spectrum is typically obtained using radiation enviroment models such as Ae8/Ap8 and Ae9/Ap9

    Returns:
        2d ndarray of the electron particle energy (column 0) vs ln(fluence) (column 1)

    """

    log_electron_integral_spectrum = np.copy(electron_integral_particle_spectrum)
    log_electron_integral_spectrum[:, 0] = electron_integral_particle_spectrum[:, 0]
    log_electron_integral_spectrum[:, 1] = np.log(electron_integral_particle_spectrum[:, 1])
    return log_electron_integral_spectrum


def get_loglog_proton_integral_spectrum(proton_integral_particle_spectrum):
    """
    Generates the ln(proton particle energy) vs ln(fluence) spectrum of the integral particle spectrum in accordance with the JPL EQGAFLUX FORTRAN program as described in the GaAs Solar Cell Radiation Handbook

    Args:
        proton_integral_particle_spectrum (ndarray) : 2d array of the integral proton particle spectrum where column 0 is the proton particle energy and column 1 is the integral fluence.  The integral particle spectrum is typically obtained using radiation enviroment models such as Ae8/Ap8 and Ae9/Ap9

    Returns:
        2d ndarray of the ln(proton particle energy) (column 0) vs ln(fluence) (column 1)

    """

    loglog_proton_integral_spectrum = np.copy(proton_integral_particle_spectrum)
    loglog_proton_integral_spectrum[:, 0] = np.log(proton_integral_particle_spectrum[:, 0])
    loglog_proton_integral_spectrum[:, 1] = np.log(proton_integral_particle_spectrum[:, 1])
    return loglog_proton_integral_spectrum


def get_loglog_rdc(relative_damage_coefficients):
    """
    Generates the ln(particle energy) vs ln(relative damage coefficient) spectrum for a give relative damage coefficient spectrum in accordance with the JPL EQGAFLUX FORTRAN program as described in the GaAs Solar Cell Radiation Handbook

    Args:
        relative_damage_coefficients (ndarray) : 2d array where column 0 contains the particle energies and column 1 contains the relative damage coefficients

    Returns:
        2d ndarray of the ln(particle energy) (column 0) vs ln(relative damage coefficient) (column 1)

    """
    loglog_rdc = np.copy(relative_damage_coefficients)
    loglog_rdc[:, 0] = np.log(relative_damage_coefficients[:, 0])
    loglog_rdc[:, 1] = np.log(relative_damage_coefficients[:, 1])
    return loglog_rdc


def compute_equivalent_fluence_for_electrons(relative_damage_coefficients_electrons,
                                             electron_integral_particle_spectrum, nstep=2):
    """
    Calculates the total 1 MeV electron fluence for a given integral electron particle spectrum and relative damage coefficients.  The function is almost a line for line translation of the JPL EQGAFLUX FORTRAN program in the GaAs Solar Cell Radiation Handbook

    Args:
        relative_damage_coefficients_electrons (ndarray) : 2d array where column 0 contains the particle energies and column 1 contains the relative damage coefficients for electrons
        electron_integral_particle_spectrum (ndarray) : 2d array of the integral electron particle spectrum where column 0 is the electron particle energy and column 1 is the integral fluence.  The integral particle spectrum is typically obtained using radiation enviroment models such as Ae8/Ap8 and Ae9/Ap9
        nstep (int) : The integration fineness and is typically 2

    Returns:
        Total 1 MeV electron fluence

    """
    interpolate_rdc = sp.interpolate.interp1d(np.log(relative_damage_coefficients_electrons[:, 0]),
                                              np.log(relative_damage_coefficients_electrons[:, 1]),
                                              fill_value='extrapolate')
    interpolate_log_electron_spectrum = sp.interpolate.interp1d(electron_integral_particle_spectrum[:, 0],
                                                                np.log(electron_integral_particle_spectrum[:, 1]))

    EQV1 = 0
    D = 0
    EMLM = get_loglog_rdc(relative_damage_coefficients_electrons)
    for i in range(len(EMLM[:, 0]) - 1):
        diff = EMLM[i + 1, 0] - EMLM[i, 0]
        delta = diff / nstep
        del2 = delta / 2
        for j in range(nstep):
            SPEC1 = (EMLM[i, 0] + delta * j)
            SPEC2 = SPEC1 + delta
            DSPEC = SPEC1 + del2
            EK = np.exp(SPEC1)
            EK1 = np.exp(SPEC2)

            if (EK1 > electron_integral_particle_spectrum[-1, 0]) or (EK < electron_integral_particle_spectrum[0, 0]):
                phi1 = 0
                phi2 = 0
                DPmax = 0
            else:
                phi1 = np.nan_to_num(np.exp(interpolate_log_electron_spectrum(EK)))
                phi2 = np.nan_to_num(np.exp(interpolate_log_electron_spectrum(EK1)))
                D = np.nan_to_num(interpolate_rdc(DSPEC))
                if D == 0:
                    D = 0
                else:
                    D = np.nan_to_num(np.exp(D))
            DPhi = phi1 - phi2
            PROD1 = DPhi * D
            EQV1 = EQV1 + PROD1
    return EQV1


def compute_equivalent_fluence_for_protons(relative_damage_coefficients_protons, proton_integral_particle_spectrum,
                                           nstep=2):
    """
    Calculates the total 1 MeV electron fluence for a given integral electron particle spectrum and relative damage coefficients.  The function is almost a line for line translation of the JPL EQGAFLUX FORTRAN program in the GaAs Solar Cell Radiation Handbook

    Args:
        relative_damage_coefficients_protons (ndarray) : 2d array where column 0 contains the particle energies and column 1 contains the relative damage coefficients for protons
        proton_integral_particle_spectrum (ndarray): 2d array of the integral proton particle spectrum where column 0 is the proton particle energy and column 1 is the integral fluence.  The integral particle spectrum is typically obtained using radiation enviroment models such as Ae8/Ap8 and Ae9/Ap9

        nstep (int) : The integration fineness and is typically 2

    Returns:
        Total 1 MeV electron fluence

    """
    interpolate_rdc = sp.interpolate.interp1d(np.log(relative_damage_coefficients_protons[:, 0]),
                                              np.log(relative_damage_coefficients_protons[:, 1]),
                                              fill_value='extrapolate')
    interpolate_loglog_proton_spectrum = sp.interpolate.interp1d(np.log(proton_integral_particle_spectrum[:, 0]),
                                                                 np.log(proton_integral_particle_spectrum[:, 1]))
    loglog_proton_particle_spectrum = get_loglog_proton_integral_spectrum(proton_integral_particle_spectrum)

    PMLN = get_loglog_rdc(relative_damage_coefficients_protons)
    EQV10P = 0
    spec1 = []
    dPhi = []
    x_new = []
    rdcMEAN = []
    dspec = []
    for i in range(len(PMLN[:, 0]) - 1):
        diff = PMLN[i + 1, 0] - PMLN[i, 0]
        delta = diff / nstep
        del2 = delta / 2
        for j in range(nstep):
            SPEC1 = (PMLN[i, 0] + delta * j)
            SPEC2 = SPEC1 + delta
            DSPEC = SPEC1 + del2
            x_new.append(SPEC1)
            if (SPEC2 > loglog_proton_particle_spectrum[-1, 0]) or (SPEC1 < loglog_proton_particle_spectrum[0, 0]):
                phi1 = 0
                phi2 = 0
                DPmax = 0
            else:

                dspec.append(DSPEC)
                phi1 = np.nan_to_num(np.exp(interpolate_loglog_proton_spectrum(SPEC1)))
                phi2 = np.nan_to_num(np.exp(interpolate_loglog_proton_spectrum(SPEC2)))
                DPmax = np.nan_to_num(interpolate_rdc(DSPEC))
                if DPmax == 0:
                    DPmax = 0
                else:
                    DPmax = np.nan_to_num(np.exp(DPmax))
            rdcMEAN.append(DPmax)
            spec1.append(SPEC1)
            DPhi = phi1 - phi2
            dPhi.append(DPhi)
            PROD1 = DPhi * DPmax
            EQV10P = EQV10P + PROD1
    return EQV10P
