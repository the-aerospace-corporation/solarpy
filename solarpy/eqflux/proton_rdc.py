import numpy as np
import scipy as sp
from scipy import interpolate
from scipy import integrate
import matplotlib.pyplot as plt
import solarpy.eqflux.electron_rdc
from solarpy import degradation_equations as d_eq
from solarpy.analyze_goodness_of_fit import analyze_goodness_of_fit
import pandas as pd
from profilehooks import profile


class proton_rdc_aero(object):
    def __init__(self, qual_data, critical_remaining_factor=None, energy_to_normalize_rdc=None,
                 density_of_shield=None, shield_thickness_cm=None, shield_range_table=None, solar_cell_range_table=None,
                 solar_cell_displacement_table=None, solar_cell_thickness_cm=None, fit_type=None, fit_parameters=None):
        if critical_remaining_factor == None:
            self.critical_remaining_factor = 0.8

        else:
            self.critical_remaining_factor = critical_remaining_factor

        if energy_to_normalize_rdc == None:
            self.energy_to_normalize_rdc = 10

        else:
            self.energy_to_normalize_rdc = energy_to_normalize_rdc

        self.qual_data = qual_data
        self.particle_energies = np.unique(qual_data[:,0])
        self.model_fit = []
        self.shield_thickness_cm = shield_thickness_cm
        self.density_of_shield = density_of_shield
        self.mass_thickness_of_shield = density_of_shield * shield_thickness_cm
        self.solar_cell_range_table = solar_cell_range_table
        self.solar_cell_displacement_table = solar_cell_displacement_table
        self.solar_cell_thickness_cm = solar_cell_thickness_cm
        self.shield_range_table = shield_range_table
        self.fit_type = fit_type
        self.fit_parameters = fit_parameters
        self.energy_vs_rdc = []
        self.energy_vs_rdc_extrapolated = []
        self.omnidirectional_shielded_rdc = []
        self.unidirectional_shielded_rdc = []

        self.update_rdcs(fit_type=fit_type, fit_parameters=fit_parameters)

    def get_energy_vs_rdc(self, fit_parameters=None, fit_type=None):
        if fit_type == 'Single':
            energy_vs_rdc = solarpy.eqflux.electron_rdc.get_energy_vs_rdc(self.qual_data, self.critical_remaining_factor,
                                                                          self.energy_to_normalize_rdc, fit_parameters=fit_parameters, return_model='y')

        elif fit_type == 'Double':
            energy_vs_rdc = solarpy.eqflux.electron_rdc.get_energy_vs_rdc_double_degradation(self.qual_data,
                                                                                             self.critical_remaining_factor,
                                                                                             self.energy_to_normalize_rdc,
                                                                                             fit_parameters=fit_parameters,
                                                                                             return_model='y')

        else:
            energy_vs_rdc = solarpy.eqflux.electron_rdc.get_energy_vs_rdc(self.qual_data, self.critical_remaining_factor,
                                                                          self.energy_to_normalize_rdc, fit_parameters=fit_parameters, return_model='y')
        self.energy_vs_rdc = energy_vs_rdc[0]
        self.model_fit = energy_vs_rdc[1]
        return self.energy_vs_rdc

    def get_energy_vs_rdc_extrapolated(self, minimum_particle_energy=1e-2, maximum_particle_energy=1e2, minimum_energy_rdc=None, maximum_energy_rdc = None,
                                       number_of_points=3000):
        # if maximum_rdc_energy is None:
        #     maximum_rdc_energy = ([50, 0.61],[100, 0.56]) #proton rdcs for GaAs Pmax-Voc at 50MeV and 100 MeV normalized for 10 MeV protons from GaAs Rad Handbook
        # if minimum_energy_rdc is None:
        #     minimum_energy_rdc = [0.02, 0.1] #proton rdcs for GaAs Pmax-Voc at 0.02 MeV normalized for 10 MeV protons from GaAs Rad Handbook
        # self.energy_vs_rdc_extrapolated = rdc.extrapolate_RDC(self.energy_vs_rdc, minimum_particle_energy,
        #                                                       maximum_particle_energy, minimum_energy_rdc, maximum_rdc_energy, number_of_points)
        self.energy_vs_rdc_extrapolated = solarpy.eqflux.electron_rdc.extrapolate_RDC_loglog(self.energy_vs_rdc, minimum_particle_energy,
                                                                                             maximum_particle_energy, minimum_energy_rdc, maximum_energy_rdc, number_of_points)

        return self.energy_vs_rdc_extrapolated

    def get_omnidirectional_shielded_rdc(self):
        omnidirectional_shielded_rdc = get_omnidirectional_and_shielded_proton_rdc(
            self.energy_vs_rdc_extrapolated, self.solar_cell_range_table, self.shield_range_table,
            self.solar_cell_displacement_table, self.solar_cell_thickness_cm, self.shield_thickness_cm)
        self.omnidirectional_shielded_rdc = omnidirectional_shielded_rdc
        return self.omnidirectional_shielded_rdc

    def get_unidirectional_shielded_rdc(self, incident_angle=0):
        omnidirectional_shielded_rdc = get_unidirectional_and_shielded_proton_rdc(
            self.energy_vs_rdc_extrapolated, self.solar_cell_range_table, self.shield_range_table,
            self.solar_cell_displacement_table, self.solar_cell_thickness_cm, self.shield_thickness_cm, incident_angle=incident_angle)
        self.omnidirectional_shielded_rdc = omnidirectional_shielded_rdc
        return self.omnidirectional_shielded_rdc


    def plot_rdc(self, rdc=None, color='blue', label=None):
        if label is None:
            label = rdc
        if rdc == 'unidirectional':
            ax = plt.loglog(self.energy_vs_rdc[:, 0], self.energy_vs_rdc[:, 1], '-o', color=color, label = label)

        elif rdc == 'extrapolated':
            ax = plt.loglog(self.energy_vs_rdc_extrapolated[:, 0], self.energy_vs_rdc_extrapolated[:, 1], '-', color=color, label = label)

        elif rdc == 'shielded':
            ax = plt.loglog(self.omnidirectional_shielded_rdc[:, 0], self.omnidirectional_shielded_rdc[:, 1], '-', color=color, label = label +' '+ str(self.shield_thickness_cm/0.00254)+'mils')
            # ax = plt.loglog(self.omnidirectional_shielded_rdc[:, 0], self.omnidirectional_shielded_rdc[:, 1], '-', color=color)

        else:
            ax = plt.loglog(self.energy_vs_rdc[:, 0], self.energy_vs_rdc[:, 1], '-o', color=color, label = 'unidirectional')

        return ax

    def plot_omnidirectional_shielded_rdc(self):
        ax = plt.loglog(self.omnidirectional_shielded_rdc[:, 0], self.omnidirectional_shielded_rdc[:, 1])
        return ax

    def update_rdcs(self, fit_type=None, fit_parameters=None, minimum_particle_energy=1e-2, maximum_particle_energy=1e2, minimum_energy_rdc=None, maximum_energy_rdc = None, number_of_points=3000): #minimum energy_rdc from GaAs Voc
        if fit_type == None:
            fit_type = self.fit_type

        if fit_parameters == None:
            fit_parameters = self.fit_parameters

        self.get_energy_vs_rdc(fit_type=fit_type, fit_parameters=fit_parameters)

        self.get_energy_vs_rdc_extrapolated(minimum_particle_energy=minimum_particle_energy,
                                            maximum_particle_energy=maximum_particle_energy,
                                            minimum_energy_rdc=minimum_energy_rdc,
                                            maximum_energy_rdc=maximum_energy_rdc,
                                            number_of_points=number_of_points)

        self.get_omnidirectional_shielded_rdc()

    def plot_to_check_fits(self, qual_data=None, fit_type=None, fit_parameters=None):
        if qual_data == None:
            qual_data=self.qual_data

        if fit_type == None:
            fit_type = self.fit_type

        if fit_parameters == None:
            fit_parameters = self.fit_parameters
        d_eq.plot_fit_checks(qual_data, fit_type, fit_parameters)

    def save_rdcs(self,file_name):
        self._save_single_rdc(file_name+'proton_rdc.txt', self.energy_vs_rdc)
        self._save_single_rdc(file_name+'proton_rdc_extrapolated.txt', self.energy_vs_rdc_extrapolated)
        shield_in_mils = self.shield_thickness_cm / 0.00254
        self._save_single_rdc(file_name+'proton_rdc_shielded_' + str(shield_in_mils) + '_mils.txt', self.omnidirectional_shielded_rdc, shielding_in_mils=shield_in_mils)
        # np.savetxt(file_name+'proton_rdc.txt', self.energy_vs_rdc, fmt = '%.6e', delimiter='\t', header='Energy (MeV)\tRDC\tRemaining Factor=\t'+str(self.critical_remaining_factor))
        # np.savetxt(file_name+'proton_rdc_extrapolated.txt', self.energy_vs_rdc_extrapolated, fmt = '%.6e', delimiter='\t', header='Energy (MeV)\tRDC\tRemaining Factor=\t'+str(self.critical_remaining_factor))
        # shield_in_mils = self.shield_thickness_cm/0.00254
        # np.savetxt(file_name+'proton_rdc_shielded_' + str(shield_in_mils) + '_mils.txt', self.omnidirectional_shielded_rdc, fmt = '%.6e', delimiter='\t', header='Energy (MeV)\tRDC\tRemaining Factor=\t'+str(self.critical_remaining_factor))

    def _save_single_rdc(self, file_name, energy_vs_rdc, shielding_in_mils=None):
        data_to_save = []
        remaining_factor_used_for_rdc = ['Remaining Factor =', self.critical_remaining_factor]
        data_to_save.append(remaining_factor_used_for_rdc)

        if shielding_in_mils:
            shield_used = ['Shielding =', 'CMG-JPL GaAs Handbook']
            shield_thickness_mils = ['Shield in mils =', shielding_in_mils]
            data_to_save.append(shield_used)
            data_to_save.append(shield_thickness_mils)
        data_to_save.append([])

        data_to_save.append(['Fit Coeffecients'])
        for i, model in enumerate(self.model_fit):
            if i == 0 and len(model.coefficients)<=2:
                data_to_save.append(['Energy (MeV)', 'C', 'phi', 'r2'])
                data_to_save.append([self.particle_energies[i]] + list(model.coefficients) + [model.r2])
            elif i == 0 and len(model.coefficients)>2:
                data_to_save.append(['Energy (MeV)', 'C1', 'phi_1',  'C2', 'phi_2', 'r2'])
                data_to_save.append([self.particle_energies[i]] + list(model.coefficients) + [model.r2])
            else:
                data_to_save.append([self.particle_energies[i]]+list(model.coefficients)+[model.r2])
        data_to_save.append([])

        data_to_save.append(['Energy (MeV)', 'RDC'])
        for row in energy_vs_rdc:
            data_to_save.append(row)

        df = pd.DataFrame(data_to_save)
        df.to_csv(file_name, sep='\t', index=False, float_format='%.8e', header=False)

def get_omnidirectional_and_shielded_proton_rdc(energy_vs_rdc, solar_cell_range_table, shield_range_table,
                                                solar_cell_displacement_table, solar_cell_thickness_cm,
                                                shield_thickness_cm):
    """
    Calculates the omnidirectional relative damage coefficient curve through shielding for protons following the procedures in the GaAs Solar Cell Radiation Handbook.  A key difference between this calculation and that for electrons is that protons can stop in the solar cell where they deposit most of their energy.  This end of track damage at the Bragg peak can lead to greater changes in the solar cell performance parameters.  The JPL method accounts for end of track damage by relating the the stopping distance of different energy particles and angles to that of the relative damage coefficient of the normally incident proton.  Also the relative damage coefficient is weighted by the ratio total displacements created for particle at an angle and its equivalent normal incidence energy at the same stopping depth.

    Args:
        energy_vs_rdc (ndarray) : 2d numpy array with column 0 being the particle energy and column 1 being the relative damage coefficient
        solar_cell_range_table (ndarray) : 2d numpy array of stopping range for protons in the solar cell material of interest.  Ideally this range table would be for the exact solar cell, but the range table for GaAs is typically used.  Column 1 is the proton energy and column 2 is the range in cm
        shield_range_table (ndarray) : 2d numpy array where column 1 is the proton energy and column 2 is the range in cm through the shielding material. For solar cells the shield is typically some type of glass.  Stopping and range table for a shield
        solar_cell_displacement_table (ndarray) : 2d numpy array where column is proton energy and column 2 is the number of displacement created by the proton in the solar cell material.  The table is generated using TRIM 2018.  By using TRIM the number of displacements can be determined for each particle energy.  In the GaAs Solar Cell handbook the displacements are calculated using a the Kinche-Pease model
        solar_cell_thickness_cm (float) : Thickness of the solar cell in cm
        shield_thickness_cm (float) : Thickness of the shield in cm


    Returns (ndarray):
        2d numpy array where column 1 is particle energy and column 2 is the relative damage coefficient

    """
    if solar_cell_range_table[0, 0] != 0:  # need 0 energy to be 0 range to get zeros for particles of low energy
        solar_cell_range_table = np.vstack(([0, 0], solar_cell_range_table))

    # if energy_vs_rdc[0, 0] != 0:  # need 0 energy get zeros for particles of low energy
    #     energy_vs_rdc = np.vstack(([0, 0], energy_vs_rdc))

    #Make sure RDC is in range of Glass and solar cell range tables
    energy_vs_rdc = energy_vs_rdc[energy_vs_rdc[:,0]<=np.max(shield_range_table[:,0])]
    energy_vs_rdc = energy_vs_rdc[energy_vs_rdc[:,0]>=np.min(shield_range_table[:,0])]
    energy_vs_rdc = energy_vs_rdc[energy_vs_rdc[:,0]<=np.max(solar_cell_range_table[:,0])]
    energy_vs_rdc = energy_vs_rdc[energy_vs_rdc[:,0]>=np.min(solar_cell_range_table[:,0])]

    anglesToIntegrateOver = np.linspace(0.0, np.pi / 2, 5000)
    newRDC = []
    intpRangeSolarCell = sp.interpolate.interp1d(solar_cell_range_table[:, 0], solar_cell_range_table[:, 1])
    for energy in energy_vs_rdc[:, 0]:
        energyAfterGlass = energyAfterGlass_single(energy, anglesToIntegrateOver, shield_thickness_cm,
                                                   shield_range_table)
        rangeInSolarCell = intpRangeSolarCell(energyAfterGlass) * np.cos(anglesToIntegrateOver)
        indicesFullPenetratingProtons = np.where(rangeInSolarCell > solar_cell_thickness_cm)
        indicesStoppedProtons = np.where(rangeInSolarCell <= solar_cell_thickness_cm)

        rdcFullPenetratingProtons = _rdcFullyPentrating(energyAfterGlass[indicesFullPenetratingProtons],
                                                        anglesToIntegrateOver[indicesFullPenetratingProtons],
                                                        energy_vs_rdc)
        rdcStoppedProtons = _rdcStoppedProtons(energyAfterGlass[indicesStoppedProtons],
                                               anglesToIntegrateOver[indicesStoppedProtons], energy_vs_rdc,
                                               solar_cell_range_table, solar_cell_displacement_table)
        if np.isnan(rdcStoppedProtons):
            rdcStoppedProtons = 0
        newRDC.append(rdcFullPenetratingProtons + rdcStoppedProtons)

    omniRDCvsEnergy = np.copy(energy_vs_rdc)
    omniRDCvsEnergy[:, 1] = np.array(newRDC)
    return omniRDCvsEnergy

def get_unidirectional_and_shielded_proton_rdc(energy_vs_rdc, solar_cell_range_table, shield_range_table,
                                               solar_cell_displacement_table, solar_cell_thickness_cm,
                                               shield_thickness_cm, incident_angle=0):
    """
    Calculates the unidirectional relative damage coefficient curve through shielding for protons following the procedures in the GaAs Solar Cell Radiation Handbook.  A key difference between this calculation and that for electrons is that protons can stop in the solar cell where they deposit most of their energy.  This end of track damage at the Bragg peak can lead to greater changes in the solar cell performance parameters.  The JPL method accounts for end of track damage by relating the the stopping distance of different energy particles and angles to that of the relative damage coefficient of the normally incident proton.  Also the relative damage coefficient is weighted by the ratio total displacements created for particle at an angle and its equivalent normal incidence energy at the same stopping depth.

    Args:
        energy_vs_rdc (ndarray) : 2d numpy array with column 0 being the particle energy and column 1 being the relative damage coefficient
        solar_cell_range_table (ndarray) : 2d numpy array of stopping range for protons in the solar cell material of interest.  Ideally this range table would be for the exact solar cell, but the range table for GaAs is typically used.  Column 1 is the proton energy and column 2 is the range in cm
        shield_range_table (ndarray) : 2d numpy array where column 1 is the proton energy and column 2 is the range in cm through the shielding material. For solar cells the shield is typically some type of glass.  Stopping and range table for a shield
        solar_cell_displacement_table (ndarray) : 2d numpy array where column is proton energy and column 2 is the number of displacement created by the proton in the solar cell material.  The table is generated using TRIM 2018.  By using TRIM the number of displacements can be determined for each particle energy.  In the GaAs Solar Cell handbook the displacements are calculated using a the Kinche-Pease model
        solar_cell_thickness_cm (float) : Thickness of the solar cell in cm
        shield_thickness_cm (float) : Thickness of the shield in cm
        incident_angle (float or ndarray): angle radiation is entering the solar cell in degrees


    Returns (ndarray):
        2d numpy array where column 1 is particle energy and column 2 is the relative damage coefficient

    """
    if solar_cell_range_table[0, 0] != 0:  # need 0 energy to be 0 range to get zeros for particles of low energy
        solar_cell_range_table = np.vstack(([0, 0], solar_cell_range_table))

    # if energy_vs_rdc[0, 0] != 0:  # need 0 energy get zeros for particles of low energy
    #     energy_vs_rdc = np.vstack(([0, 0], energy_vs_rdc))

    #Make sure RDC is in range of Glass and solar cell range tables
    energy_vs_rdc = energy_vs_rdc[energy_vs_rdc[:,0]<=np.max(shield_range_table[:,0])]
    energy_vs_rdc = energy_vs_rdc[energy_vs_rdc[:,0]>=np.min(shield_range_table[:,0])]
    energy_vs_rdc = energy_vs_rdc[energy_vs_rdc[:,0]<=np.max(solar_cell_range_table[:,0])]
    energy_vs_rdc = energy_vs_rdc[energy_vs_rdc[:,0]>=np.min(solar_cell_range_table[:,0])]

    anglesToIntegrateOver = incident_angle
    if isinstance(anglesToIntegrateOver, np.ndarray) is not True:
        anglesToIntegrateOver = np.deg2rad(np.array([anglesToIntegrateOver]))
    else:
        anglesToIntegrateOver = np.deg2rad(anglesToIntegrateOver)

    newRDC = []
    intpRangeSolarCell = sp.interpolate.interp1d(solar_cell_range_table[:, 0], solar_cell_range_table[:, 1])
    for energy in energy_vs_rdc[:, 0]:
        energyAfterGlass = energyAfterGlass_single(energy, anglesToIntegrateOver, shield_thickness_cm,
                                                   shield_range_table)
        rangeInSolarCell = intpRangeSolarCell(energyAfterGlass) * np.cos(anglesToIntegrateOver)
        indicesFullPenetratingProtons = np.where(rangeInSolarCell > solar_cell_thickness_cm)
        indicesStoppedProtons = np.where(rangeInSolarCell <= solar_cell_thickness_cm)
        rdcFullPenetratingProtons = _shieldedRDC_Single(energyAfterGlass[indicesFullPenetratingProtons], energy_vs_rdc)
        rdcStoppedProtons = _rdcStoppedProtons(energyAfterGlass[indicesStoppedProtons],
                                               anglesToIntegrateOver[indicesStoppedProtons], energy_vs_rdc,
                                               solar_cell_range_table, solar_cell_displacement_table, hemisphere=False)
        if np.isnan(rdcStoppedProtons)or not rdcStoppedProtons:
            rdcStoppedProtons = 0
        if np.isnan(rdcFullPenetratingProtons)or not rdcFullPenetratingProtons:
            rdcFullPenetratingProtons = 0
        newRDC.append(rdcFullPenetratingProtons + rdcStoppedProtons)

    omniRDCvsEnergy = np.copy(energy_vs_rdc)
    omniRDCvsEnergy[:, 1] = np.array(newRDC, dtype=object).flatten()
    return omniRDCvsEnergy


def _rdcFullyPentrating(particle_energy, angles, energy_vs_rdc):
    shieldedRDCs = _shieldedRDC_Single(particle_energy, energy_vs_rdc)
    if not shieldedRDCs.any():
        shieldedRDCs = 0
    else:
        shieldedRDCs = shieldedRDCs
    y = shieldedRDCs * 2 * np.pi * np.sin(angles)
    rdcFullyPenetrating = (1 / (4 * np.pi)) * (sp.integrate.trapz(y, angles))
    return rdcFullyPenetrating

def _rdcStoppedProtons(particle_energy, angles, energy_vs_rdc, solar_cell_range_table, solar_cell_displacement_table, hemisphere=True):
    if solar_cell_range_table[0, 0] != 0:  # need 0 energy to be 0 range to get zeros for particles of low energy
        solar_cell_range_table = np.vstack(([0, 0], solar_cell_range_table))

    if solar_cell_displacement_table[0, 0] != 0:  # need 0 energy to be 0 range to get zeros for particles of low energy
        solar_cell_displacement_table = np.vstack(([0, 0], solar_cell_displacement_table))

    if energy_vs_rdc[0, 0] != 0:  # need 0 energy to be 0 range to get zeros for particles of low energy
        energy_vs_rdc = np.vstack(([0, 0], energy_vs_rdc))

    # intpRDC = sp.interpolate.interp1d(energy_vs_rdc[:, 0], energy_vs_rdc[:, 1])
    # intpRangeSolarCell = sp.interpolate.interp1d(solar_cell_range_table[:, 0], solar_cell_range_table[:, 1])
    # intpEnergySolarCell = sp.interpolate.interp1d(solar_cell_range_table[:, 1], solar_cell_range_table[:, 0])
    # intpDisplacements = sp.interpolate.interp1d(solar_cell_displacement_table[:, 0], solar_cell_displacement_table[:, 1])

    # totalDisplacementsFromParticlesAfterGlass = intpDisplacements(particle_energy)
    totalDisplacementsFromParticlesAfterGlass = np.interp(particle_energy, solar_cell_displacement_table[:, 0], solar_cell_displacement_table[:, 1])

    # normalIncidenceParticleEnergies = intpEnergySolarCell(intpRangeSolarCell(particle_energy) * np.cos(angles))
    interpolated_range = np.interp(particle_energy, solar_cell_range_table[:, 0], solar_cell_range_table[:, 1]) *np.cos(angles)
    normalIncidenceParticleEnergies = np.interp(interpolated_range, solar_cell_range_table[:, 1], solar_cell_range_table[:, 0])
    # totalDisplacementsFromNormalIncidenceParticles = intpDisplacements(normalIncidenceParticleEnergies)
    totalDisplacementsFromNormalIncidenceParticles = np.interp(normalIncidenceParticleEnergies, solar_cell_displacement_table[:, 0], solar_cell_displacement_table[:, 1])

    # rdcOfNormalIncidenceParticles = intpRDC(normalIncidenceParticleEnergies)
    rdcOfNormalIncidenceParticles = np.interp(normalIncidenceParticleEnergies, energy_vs_rdc[:, 0], energy_vs_rdc[:, 1])
    rdcStoppedProtons = rdcOfNormalIncidenceParticles * (totalDisplacementsFromParticlesAfterGlass / totalDisplacementsFromNormalIncidenceParticles)
    rdcStoppedProtons = np.nan_to_num(rdcStoppedProtons)
    if hemisphere:
        y = rdcStoppedProtons * 2 * np.pi * np.sin(angles) * np.cos(angles)
        rdcStoppedProtons = (1 / (4 * np.pi)) * (sp.integrate.trapz(y, angles))

    else:
        rdcStoppedProtons = rdcStoppedProtons * np.cos(angles) # have to include projected area for stopped protons
    return rdcStoppedProtons


def energyAfterGlass_single(unidirectional_particle_energy, entry_angles, shield_thickness_cm, shield_range_table):
    if shield_range_table[0, 0] != 0:  # need 0 energy to be 0 range to get zeros for particles of low energy
        shield_range_table = np.vstack(([0, 0], shield_range_table))

    # interpolate_Energy = sp.interpolate.interp1d(shield_range_table[:, 1], shield_range_table[:, 0])
    shieldedRange = get_range_over_entry_angles(unidirectional_particle_energy, entry_angles, shield_thickness_cm,
                                                shield_range_table)
    # shieldedParticleEnergies = interpolate_Energy(shieldedRange)
    shieldedParticleEnergies = np.interp(shieldedRange, shield_range_table[:, 1], shield_range_table[:, 0])
    return np.array(shieldedParticleEnergies)


def _shieldedRDC_Single(particle_energy, energy_vs_rdc):
    intpRDC = sp.interpolate.interp1d(energy_vs_rdc[:, 0], energy_vs_rdc[:, 1], bounds_error=False, fill_value=0)
    if particle_energy[particle_energy < 0].any():
        particle_energy *= (particle_energy > 0)
    shieldedRDCs = intpRDC(particle_energy)


    return shieldedRDCs


def get_range_over_entry_angles(unidirectional_particle_energy, entry_angles, shield_thickness_cm,
                                shield_range_table):
    # interpolate_Range = sp.interpolate.interp1d(shield_range_table[:, 0], shield_range_table[:, 1])
    # range_over_entry_angles = interpolate_Range(unidirectional_particle_energy) - (shield_thickness_cm / np.cos(entry_angles))
    range_over_entry_angles = np.interp(unidirectional_particle_energy, shield_range_table[:, 0], shield_range_table[:, 1]) - (shield_thickness_cm / np.cos(entry_angles))
    range_over_entry_angles *= (range_over_entry_angles > 0)  # Gets rid of negative numbers and make them 0....because there is not such think as a negative range
    return np.array(range_over_entry_angles)


def _totalNumberOfDisplacements(particle_energy, solar_cell_displacement_table):
    intpDisplacements = sp.interpolate.interp1d(solar_cell_displacement_table[:, 0],
                                                solar_cell_displacement_table[:, 1])
    totalDisplacements = intpDisplacements(particle_energy)
    return totalDisplacements
