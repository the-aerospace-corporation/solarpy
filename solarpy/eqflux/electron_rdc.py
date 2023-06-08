import numpy as np
import scipy as sp
from scipy import interpolate
from scipy import integrate
from scipy import stats
from solarpy import degradation_equations as d_eq
from solarpy.analyze_goodness_of_fit import analyze_goodness_of_fit
import matplotlib.pyplot as plt
import pandas as pd

np.seterr(divide='ignore', invalid='ignore')

class electron_rdc_aero(object):
    def __init__(self, qual_data, critical_remaining_factor=None, energy_to_normalize_rdc=None, density_of_shield=None,
                 shield_thickness_cm=None, shield_range_table=None, fit_type=None, fit_parameters=None, normalize_rdc=False):
        if critical_remaining_factor is None:
            self.critical_remaining_factor = 0.8

        else:
            self.critical_remaining_factor = critical_remaining_factor

        if energy_to_normalize_rdc is None:
            self.energy_to_normalize_rdc = 1

        else:
            self.energy_to_normalize_rdc = energy_to_normalize_rdc

        self.qual_data = qual_data
        self.particle_energies = np.unique(qual_data[:,0])
        self.model_fit = []
        self.shield_thickness_cm = shield_thickness_cm
        self.density_of_shield = density_of_shield
        self.mass_thickness_of_shield = density_of_shield * shield_thickness_cm
        self.shield_range_table = shield_range_table
        self.fit_type = fit_type
        self.fit_parameters= fit_parameters
        self.energy_vs_rdc = []
        self.energy_vs_rdc_extrapolated = []
        self.omnidirectional_shielded_rdc = []
        self.unidirectional_shielded_rdc = []
        self.normalize_rdc = normalize_rdc

        self.update_rdcs()

    def get_energy_vs_rdc(self, fit_parameters=None, fit_type=None):

        if fit_type.lower() == 'single':
            energy_vs_rdc = get_energy_vs_rdc(self.qual_data,self.critical_remaining_factor,self.energy_to_normalize_rdc, fit_parameters=fit_parameters, return_model='y')

        elif fit_type.lower() == 'double':
            energy_vs_rdc = get_energy_vs_rdc_double_degradation(self.qual_data, self.critical_remaining_factor, self.energy_to_normalize_rdc, fit_parameters=fit_parameters, return_model='y')

        else:
            energy_vs_rdc = get_energy_vs_rdc(self.qual_data, self.critical_remaining_factor, self.energy_to_normalize_rdc, fit_parameters=fit_parameters, return_model='y')

        self.energy_vs_rdc = energy_vs_rdc[0]
        self.model_fit = energy_vs_rdc[1]
        return self.energy_vs_rdc

    def get_energy_vs_rdc_extrapolated(self, minimum_particle_energy=1e-5, maximum_particle_energy=1e2, minimum_energy_rdc=None, maximum_energy_rdc=None,
                                       number_of_points=3000):
        if minimum_energy_rdc is None:
            minimum_energy_rdc = [0.260, 0.1] #electron rdcs for GaAs at 0.260MeV normalized for 1 MeV electrons from GaAs Rad Handbook
        # self.energy_vs_rdc_extrapolated = rdc.extrapolate_RDC(self.energy_vs_rdc,
        #                                                       minimum_particle_energy, maximum_particle_energy, minimum_energy_rdc, maximum_rdc_energy,
        #                                                       number_of_points)
        self.energy_vs_rdc_extrapolated = extrapolate_RDC_loglog(self.energy_vs_rdc, minimum_particle_energy,
                                                                                             maximum_particle_energy, minimum_energy_rdc, maximum_energy_rdc, number_of_points)

        self.energy_vs_rdc_extrapolated[self.energy_vs_rdc_extrapolated[:,0]<minimum_energy_rdc[0],1]=0
        return self.energy_vs_rdc_extrapolated

    def get_omnidirectional_shielded_rdc(self):
        self.mass_thickness_of_shield = self.density_of_shield*self.shield_thickness_cm
        omnidirectional_shielded_rdc = get_omnidirectional_electron_RDC_with_shielding(self.energy_vs_rdc_extrapolated,
                                                                                       self.mass_thickness_of_shield,
                                                                                       self.shield_range_table)
        self.omnidirectional_shielded_rdc = omnidirectional_shielded_rdc
        return self.omnidirectional_shielded_rdc

    def get_unidirectional_shielded_rdc(self, incident_angle=0):
        self.mass_thickness_of_shield = self.density_of_shield * self.shield_thickness_cm
        omnidirectional_shielded_rdc = get_unidirectional_electron_RDC_with_shielding(self.energy_vs_rdc_extrapolated,
                                                                                      self.mass_thickness_of_shield,
                                                                                      self.shield_range_table,
                                                                                      incident_angle)
        self.omnidirectional_shielded_rdc = omnidirectional_shielded_rdc
        return self.omnidirectional_shielded_rdc


    def plot_rdc(self, rdc=None, color='blue'):
        print(rdc)
        if rdc == 'unidirectional':
            ax = plt.loglog(self.energy_vs_rdc[:, 0], self.energy_vs_rdc[:, 1], '-o', color=color, label = 'unidirectional')

        elif rdc == 'extrapolated':
            ax = plt.loglog(self.energy_vs_rdc_extrapolated[:, 0], self.energy_vs_rdc_extrapolated[:, 1], '-', color=color, label = 'extrapolated')

        elif rdc == 'shielded':
            ax = plt.loglog(self.omnidirectional_shielded_rdc[:, 0], self.omnidirectional_shielded_rdc[:, 1], '-', color=color, label = 'shielded '+ str(self.shield_thickness_cm/0.00254)+'mils')

        else:
            ax = plt.loglog(self.energy_vs_rdc[:, 0], self.energy_vs_rdc[:, 1], '-o', color=color, label = 'unidirectional')

        return ax

    def plot_omnidirectional_shielded_rdc(self):
        ax = plt.loglog(self.omnidirectional_shielded_rdc[:, 0], self.omnidirectional_shielded_rdc[:, 1])
        return ax

    def update_rdcs(self, fit_type=None, fit_parameters=None, minimum_particle_energy=1e-5, maximum_particle_energy=1e2, minimum_energy_rdc=None, maximum_energy_rdc=None, number_of_points=3000):
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

    def plot_to_check_fits(self, qual_data=None, fit_type=None, fit_parameters=None, normalize_rdc=None):
        if qual_data == None:
            qual_data=self.qual_data

        if fit_type == None:
            fit_type = self.fit_type

        if fit_parameters == None:
            fit_parameters = self.fit_parameters

        if normalize_rdc == None:
            normalize_rdc = self.normalize_rdc

        d_eq.plot_fit_checks(qual_data, fit_type, fit_parameters, normalize_rdc)

    def save_rdcs(self,file_name):
        self._save_single_rdc(file_name+'electron_rdc.txt', self.energy_vs_rdc)
        self._save_single_rdc(file_name+'electron_rdc_extrapolated.txt', self.energy_vs_rdc_extrapolated)
        shield_in_mils = self.shield_thickness_cm / 0.00254
        self._save_single_rdc(file_name+'electron_rdc_shielded_' + str(shield_in_mils) + '_mils.txt', self.omnidirectional_shielded_rdc, shielding_in_mils=shield_in_mils)
        # np.savetxt(file_name+'proton_rdc.txt', self.energy_vs_rdc, fmt = '%.6e', delimiter='\t', header='Energy (MeV)\tRDC\tRemaining Factor=\t'+str(self.critical_remaining_factor))
        # np.savetxt(file_name+'proton_rdc_extrapolated.txt', self.energy_vs_rdc_extrapolated, fmt = '%.6e', delimiter='\t', header='Energy (MeV)\tRDC\tRemaining Factor=\t'+str(self.critical_remaining_factor))
        # shield_in_mils = self.shield_thickness_cm/0.00254
        # np.savetxt(file_name+'proton_rdc_shielded_' + str(shield_in_mils) + '_mils.txt', self.omnidirectional_shielded_rdc, fmt = '%.6e', delimiter='\t', header='Energy (MeV)\tRDC\tRemaining Factor=\t'+str(self.critical_remaining_factor))

    def _save_single_rdc(self, file_name, energy_vs_rdc, shielding_in_mils=None):
        data_to_save = []
        remaining_factor_used_for_rdc = ['Remaining Factor=', self.critical_remaining_factor]
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


def get_energy_vs_rdc(qual_data, critical_remaining_factor, energyToNormalizeRDC=1, fit_parameters=None, return_model=None):
    """
    Calculates the unidirectional RDC curve (energy vs RDC) given the radiation qual data, critical remaining factor, and the partcle energy you wish to relate the damage of each particel too. Uses the single degradation equation.  This is needed to calculate the fluence at which each the critical remaining factor is for each of the particle energy data sets

    Args:
        qual_data: 2d numpy array with three columns where column 0 is particle energy, column 1 is fluence, and column 2 is remaining factor
        critical_remaining_factor: The remaining factor at which to relate all fluences
        energyToNormalizeRDC: The particle energy to normalize the relative damage coefficients
        fit_parameters: 1d numpy array include 2 parameters, one for the fitting constant C and the other for phi_x.
    
    Returns:
        2d numpy array of the energy in column 0 and relative damage coefficient in column 1

    """
    particleEnergies = np.unique(qual_data[:, 0])
    indexOfEnergyToNormalizeRDC = np.where(particleEnergies == energyToNormalizeRDC)
    qualDataGroupedByParticleEnergy = group_qual_data_by_particle_energy(qual_data)
    critical_fluences = []
    model_list = []

    if any(isinstance(param, list) for param in fit_parameters):
        multi_param = True
    else:
        multi_param = False

    for i, qualData in enumerate(qualDataGroupedByParticleEnergy):
        fluence_vs_remaining_factor = qualData[:, [1, 2]]
        model = d_eq.singleDegradationEquation(fluence_vs_remaining_factor)
        if multi_param:
            fit_param = fit_parameters[i]
        else:
            fit_param = fit_parameters

        model.fit(fit_param)
        model_list.append(model)
        # print model.r2
        # fit_coefficients = model.coefficients
        # fit_coefficients = d_eq.fit_degradation_equation(fluence_vs_remaining_factor, fit_parameters)
        critical_fluence = model.getFluence(critical_remaining_factor)
        critical_fluences.append(critical_fluence)

    critical_fluences = np.array(critical_fluences)
    energyVSrdc = np.zeros((len(critical_fluences), 2))
    energyVSrdc[:, 0] = np.unique(particleEnergies)
    energyVSrdc[:, 1] = critical_fluences[indexOfEnergyToNormalizeRDC] / critical_fluences

    if return_model in ['y', 'Y']:
        return [energyVSrdc, model_list]

    else:
        return energyVSrdc


def get_energy_vs_rdc_double_degradation(qual_data, critical_remaining_factor, energyToNormalizeRDC=1, fit_parameters=None, return_model=None):
    """
    Calculates the unidirectional RDC curve (energy vs RDC) given the radiation qual data, critical remaining factor, and the partcle energy you wish to relate the damage of each particel too. Uses double degradatin equation. This is needed to calculate the fluence at which each the critical remaining factor is for each of the particle energy data sets

    Args:
        qual_data: 2d numpy array with three columns where column 0 is particle energy, column 1 is fluence, and column 2 is remaining factor
        critical_remaining_factor: The remaining factor at which to relate all fluences
        energyToNormalizeRDC: The particle energy to normalize the relative damage coefficients
        fit_parameters:  1d numpy array of list of parameters includes 4 parameters, C_1 and phi_x_1 for the first degradation equation and C_2 and phi_x_2 for the second degradation equation

    Returns:
        2d numpy array of the energy in column 0 and relative damage coefficient in column 1

    """
    if any(isinstance(param, list) for param in fit_parameters):
        multi_param = True
    else:
        multi_param = False

    particleEnergies = np.unique(qual_data[:, 0])
    indexOfEnergyToNormalizeRDC = np.where(particleEnergies == energyToNormalizeRDC)
    qualDataGroupedByParticleEnergy = group_qual_data_by_particle_energy(qual_data)
    critical_fluences = []
    model_list = []
    for i, qualData in enumerate(qualDataGroupedByParticleEnergy):
        fluence_vs_remaining_factor = qualData[:, [1, 2]]
        model = d_eq.doubleDegradationEquation(fluence_vs_remaining_factor)
        if multi_param:
            fit_param = fit_parameters[i]
        else:
            fit_param = fit_parameters
        model.fit(fit_param)
        model_list.append(model)
        critical_fluence = model.getFluence(critical_remaining_factor)
        critical_fluences.append(critical_fluence)

    critical_fluences = np.array(critical_fluences)
    energyVSrdc = np.zeros((len(critical_fluences), 2))
    energyVSrdc[:, 0] = np.unique(particleEnergies)
    energyVSrdc[:, 1] = critical_fluences[indexOfEnergyToNormalizeRDC] / critical_fluences

    if return_model in ['y', 'Y']:
        return [energyVSrdc, model_list]

    else:
        return energyVSrdc


def fitQualData(qual_data, fit_type='Single Degradation', fit_parameters=None):
    # TODO: overly complicated, may need to move this to a class of "best practices" or "tests"
    qualDataGroupedByParticelEnergy = group_qual_data_by_particle_energy(qual_data)
    data = type('data', (object,), {'fits': None, 'r2': None, 'chi2': None})
    fits = []
    rSquared = []
    chiSquared = []
    for qualData in qualDataGroupedByParticelEnergy:
        fluencVSRemainingFactor = qualData[:, [1, 2]]
        if fit_type == 'Single Degradation':
            fit = d_eq.fit_degradation_equation(fluence_or_ddd_vs_remaining_factor=fluencVSRemainingFactor,
                                                parameters=fit_parameters)
            observed = d_eq.degradation_equation(fluencVSRemainingFactor[:, 0], *fit)
            r2 = analyze_goodness_of_fit(fluencVSRemainingFactor[:, 1], observed).r2
            chi2 = sp.stats.chisquare(observed, fluencVSRemainingFactor[:, 1])[0]

        elif fit_type == 'Double Degradation':
            fit = d_eq.fit_double_degradation_equation(fluence_or_ddd_vs_remaining_factor=fluencVSRemainingFactor,
                                                       parameters=fit_parameters)
            observed = d_eq.double_degradation_equation(fluencVSRemainingFactor[:, 0], *fit)
            r2 = analyze_goodness_of_fit(fluencVSRemainingFactor[:, 1], observed).r2
            chi2 = sp.stats.chisquare(observed, fluencVSRemainingFactor[:, 1])[0]

        else:
            ###Default to single degradation model###
            fit = d_eq.fit_degradation_equation(fluence_or_ddd_vs_remaining_factor=fluencVSRemainingFactor,
                                                parameters=fit_parameters)
            observed = d_eq.degradation_equation(fluencVSRemainingFactor[:, 0], *fit)
            r2 = analyze_goodness_of_fit(fluencVSRemainingFactor[:, 1], observed).r2
            chi2 = sp.stats.chisquare(observed, fluencVSRemainingFactor[:, 1])[0]

        fits.append(fit)
        rSquared.append(r2)
        chiSquared.append(chi2)
    data.fits = fits
    data.r2 = rSquared
    data.chi2 = chiSquared
    return data


def extrapolate_RDC(energy_vs_rdc, minimum_particle_energy, maximum_particle_energy, minimum_energy_rdc = None, maximum_energy_rdc = None, number_of_points=3000):
    """
    Linearly extrapolates AND interpolates the RDC data. This is a wrapper function of scipy interpolate using the 'extrapolate" keyword argument. The user can use any function they wish to interpolate or extrapolate

    Args:
        energy_vs_rdc: 2d numpy array with column 0 being the particle energy and column 1 being the relative damage coefficienct
        minimum_particle_energy: Minimum energy to start linear extrapolation of rdc data
        maximum_particle_energy: Maximum energy to linearly extrapolate rdc data.  For example the GaAs data from the GaAs Solar Cell Radiation handbook extrapolates the low energy electron data and the data from 12MeV electrons to 40MeV electrons is extrapoloated
        minimum_energy_rdc:  1d array of minimum energy and rdc value to extrapolate.  The first value is the particle energy and the second value is the rdc.  For electrons, the minimum energy used in The GaAs Solar Cell Radiation Handbook is 0.260 MeV and the minimum rdc value for that energy is <0.1 (B. Anspaugh, 'Proton and electron damage coefficients for GaAs/Ge solar cells', DOI: 10.1109/PVSC.1991.169472)
        maximum_energy_rdc: 1d array of maximum energy and rdc value to extrapolate.  For most new triple junction cells peeople use the RDCs greater than 10 MeV from the GaAs RDC data in The GaAs Solar Cell Radiation Handbook....which...I don't agree with because there is no reference to the data see (B. Anspaugh, 'Proton and electron damage coefficients for GaAs/Ge solar cells', DOI: 10.1109/PVSC.1991.169472)
        number_of_points: The number of particle energies to interpolate. refer to numpy.logspace documentation to
        understand more

    Returns:
        2d numpy array or energy in column 0 and relative damage coefficients in column 2 that is extrapolated and interpolated for user defined arguments

    """
    if minimum_energy_rdc is not None:
        energy_vs_rdc = np.vstack((minimum_energy_rdc, energy_vs_rdc))

    if maximum_energy_rdc is not None:
        energy_vs_rdc = np.vstack((energy_vs_rdc, maximum_energy_rdc))

    if energy_vs_rdc[0,0] != 0:
        energy_vs_rdc = np.vstack(([0,0], energy_vs_rdc))

    intpRDC = sp.interpolate.interp1d(energy_vs_rdc[:, 0], energy_vs_rdc[:, 1], fill_value='extrapolate')
    particle_energies_i = np.logspace(np.log10(minimum_particle_energy), np.log10(maximum_particle_energy),
                                      number_of_points)
    rdc_i = intpRDC(particle_energies_i)
    energyVSrdc_i = np.zeros((len(particle_energies_i), 2))
    energyVSrdc_i[:, 0] = particle_energies_i
    energyVSrdc_i[:, 1] = np.nan_to_num(rdc_i)
    return energyVSrdc_i


def extrapolate_RDC_loglog(energy_vs_rdc, minimum_particle_energy, maximum_particle_energy, minimum_energy_rdc = None, maximum_energy_rdc = None, number_of_points=3000):
    """
    Linearly extrapolates AND interpolates the loglog of RDC data. This is a wrapper function of scipy interpolate using the 'extrapolate" keyword argument. The user can use any function they wish to interpolate or extrapolate

    Args:
        energy_vs_rdc: 2d numpy array with column 0 being the particle energy and column 1 being the relative damage coefficienct
        minimum_particle_energy: Minimum energy to start linear extrapolation of rdc data
        maximum_particle_energy: Maximum energy to linearly extrapolate rdc data.  For example the GaAs data from the GaAs Solar Cell Radiation handbook extrapolates the low energy electron data and the data from 12MeV electrons to 40MeV electrons is extrapoloated
        minimum_energy_rdc:  1d array of minimum energy and rdc value to extrapolate.  The first value is the particle energy and the second value is the rdc.  For electrons, the minimum energy used in The GaAs Solar Cell Radiation Handbook is 0.260 MeV and the minimum rdc value for that energy is <0.1 (B. Anspaugh, 'Proton and electron damage coefficients for GaAs/Ge solar cells', DOI: 10.1109/PVSC.1991.169472).
        maximum_energy_rdc: 1d array of maximum energy and rdc value to extrapolate.  For most new triple junction cells peeople use the RDCs greater than 10 MeV from the GaAs RDC data in The GaAs Solar Cell Radiation Handbook....which...I don't agree with because there is no reference to the data see (B. Anspaugh, 'Proton and electron damage coefficients for GaAs/Ge solar cells', DOI: 10.1109/PVSC.1991.169472)
        number_of_points: The number of particle energies to interpolate. refer to numpy.logspace documentation to
        understand more

    Returns:
        2d numpy array or energy in column 0 and relative damage coefficients in column 2 that is extrapolated and interpolated for user defined arguments

    """
    if minimum_energy_rdc is not None:
        energy_vs_rdc = np.vstack((minimum_energy_rdc, energy_vs_rdc))

    if maximum_energy_rdc is not None:
        energy_vs_rdc = np.vstack((energy_vs_rdc, maximum_energy_rdc))

    # if energy_vs_rdc[0,0] != 0: # rdc curve should definitely be 0 when there is a 0 energy particle
    #     energy_vs_rdc = np.vstack(([0,0], energy_vs_rdc))

    ## log10 of Energy vs RDC
    energy_vs_rdc_loglog = np.copy(energy_vs_rdc)
    energy_vs_rdc_loglog[:,0] = np.log10(energy_vs_rdc[:,0])
    energy_vs_rdc_loglog[:,1] = np.log10(energy_vs_rdc[:,1])
    # energy_vs_rdc_loglog = energy_vs_rdc_loglog[~np.isinf(energy_vs_rdc_loglog).any(axis=1)]
    # energy_vs_rdc_loglog = np.nan_to_num(energy_vs_rdc_loglog)
    # for x in energy_vs_rdc:
    #     print(x)
    intpRDC = sp.interpolate.interp1d(energy_vs_rdc_loglog[:, 0], energy_vs_rdc_loglog[:, 1], kind='linear', fill_value='extrapolate')
    particle_energies_i = np.logspace(np.log10(minimum_particle_energy), np.log10(maximum_particle_energy),
                                      number_of_points)
    rdc_i = intpRDC(np.log10(particle_energies_i))
    energyVSrdc_i = np.zeros((len(particle_energies_i), 2))
    energyVSrdc_i[:, 0] = particle_energies_i
    energyVSrdc_i[:, 1] = np.nan_to_num(10**rdc_i)
    return energyVSrdc_i

def get_omnidirectional_electron_RDC_with_shielding(energy_vs_rdc, mass_thickness_of_shield, shield_range_table):
    """
    Calculates the omnidirectional relative damage coefficient of a solar cell after some shielding, typically solar cell cover glass.

    Args:
        energy_vs_rdc: 2d numpy array with column 0 being the particle energy and column 1 being the relative damage coefficienct
        mass_thickness_of_shield: Mass thickness of the shielding being used in units of g/cm.  Derived by multiplying the density of the shield by the thickness of the shield
        shield_range_table: Stopping power and range table for the shield

    Returns:
        2d numpy array of energy (column 0) and relative damage coefficients (column 2).  The relative damage coefficients are calculated after going through some shielded and adjusted from unidirectional test energies to omnidirectional test energies

    """
    anglesToIntegrateOver = np.linspace(0.0, np.pi / 2, 100)  # 180 degree hemisphere
    shieldedRDC = get_rdc_after_shielding(energy_vs_rdc, anglesToIntegrateOver, mass_thickness_of_shield, shield_range_table)
    newRDC = []
    for RDC in shieldedRDC:
        y = RDC * (2 * np.pi * np.sin(anglesToIntegrateOver))
        integralOfY = (1 / (4 * np.pi)) * (sp.integrate.trapz(y, anglesToIntegrateOver))
        newRDC.append(integralOfY)
    omniRDCvsEnergy = np.copy(energy_vs_rdc)
    omniRDCvsEnergy[:, 1] = np.array(np.nan_to_num(newRDC))
    return omniRDCvsEnergy

def get_unidirectional_electron_RDC_with_shielding(energy_vs_rdc, mass_thickness_of_shield, shield_range_table, incident_angle=0):
    """
    Calculates the omnidirectional relative damage coefficient of a solar cell after some shielding, typically solar cell cover glass.

    Args:
        energy_vs_rdc: 2d numpy array with column 0 being the particle energy and column 1 being the relative damage coefficienct
        mass_thickness_of_shield: Mass thickness of the shielding being used in units of g/cm.  Derived by multiplying the density of the shield by the thickness of the shield
        shield_range_table: Stopping power and range table for the shield

    Returns:
        2d numpy array of energy (column 0) and relative damage coefficients (column 2).  The relative damage coefficients are calculated after going through some shielded and adjusted from unidirectional test energies to omnidirectional test energies

    """

    anglesToIntegrateOver = incident_angle
    if isinstance(anglesToIntegrateOver, np.ndarray) is not True:
        anglesToIntegrateOver = np.array([anglesToIntegrateOver])

    shieldedRDC = get_rdc_after_shielding(energy_vs_rdc, anglesToIntegrateOver, mass_thickness_of_shield, shield_range_table)
    newRDC = np.array(shieldedRDC)[:,0]
    # for RDC in shieldedRDC:
    #     y = RDC * (2 * np.pi * np.sin(anglesToIntegrateOver))
    #     integralOfY = (1 / (4 * np.pi)) * (sp.integrate.trapz(y, anglesToIntegrateOver))
    #     newRDC.append(integralOfY)
    omniRDCvsEnergy = np.copy(energy_vs_rdc)
    omniRDCvsEnergy[:, 1] = np.array(np.nan_to_num(newRDC))
    return omniRDCvsEnergy

def get_rdc_after_shielding(energy_vs_rdc, angles, mass_thickness_of_shield, shield_range_table):
    """
    Calculates the relative damage coefficient of unidirectional particles after some shielding.  The relative damage coefficienct is calculated for all angles of a particle after shielding over a hemisphere

    Args:
        energy_vs_rdc: 2d numpy array with column 0 being the particle energy and column 1 being the relative damage coefficient
        angles: 1d numpy array of entry angles of particles. Typically over 180 degrees
        mass_thickness_of_shield: Mass thickness of the shielding being used in units of g/cm.  Derived by multiplying the density of the shield by the thickness of the shield
        shield_range_table: Stopping power and range table for the shield

    Returns:
        List of 1d numpy array of relative damage coefficients for each angle of a hemisphere after a user defined shield

    """
    intpRDC = sp.interpolate.interp1d(energy_vs_rdc[:, 0], energy_vs_rdc[:, 1], fill_value='extrapolate')
    shieldedParticleEnergies = get_shielded_particle_energies(energy_vs_rdc[:, 0],
                                                              angles, mass_thickness_of_shield, shield_range_table)
    shieldedRDCs = []
    for i, shieldedEnergy in enumerate(shieldedParticleEnergies):
        shieldedEnergy *= (shieldedEnergy>0)
        shieldedRDCs.append(intpRDC(shieldedEnergy))
    return shieldedRDCs


def get_shielded_particle_energies(particleEnergies, angles, mass_thickness_of_shield, shield_range_table):
    """
    Calculates the energy of a particle after it passes through shielding.  The energy of the particle is calculated for all angles over 180 degrees.

    Args:
        particleEnergies:
        angles:
        mass_thickness_of_shield:
        shield_range_table:

    Returns:
        List of 1d numpy array of particle energies for each angle of a hemisphere after a user defined shield

    """
    if shield_range_table[0,0] != 0: # need 0 energy to be 0 range to get zeros for particles of low energy
        shield_range_table = np.vstack(([0,0], shield_range_table))

    intpRange = sp.interpolate.interp1d(shield_range_table[:, 0], shield_range_table[:, 1])
    intpEnergy = sp.interpolate.interp1d(shield_range_table[:, 1], shield_range_table[:, 0])
    shielded_particle_energies = []
    for energy in particleEnergies:
        newRange = intpRange(energy) - (mass_thickness_of_shield / np.cos(angles))
        newRange *= (newRange > 0)
        shieldedParticle = intpEnergy(newRange)
        shielded_particle_energies.append(shieldedParticle)
    return shielded_particle_energies


def _get_indices_of_duplicates_of_particle_energy(particle_energies):
    indicesOfDuplicates = []
    energies = np.unique(particle_energies)
    for energy in energies:
        indicesOfDuplicates.append(list(np.where(particle_energies == energy)))
    return indicesOfDuplicates


def group_qual_data_by_particle_energy(qual_data):
    particleEnergies = qual_data[:, 0]
    groupedByEnergy = []
    indicesOfDuplicates = _get_indices_of_duplicates_of_particle_energy(particle_energies=particleEnergies)
    for index in indicesOfDuplicates:
        groupedByEnergy.append(qual_data[np.min(index):np.max(index) + 1])
    return groupedByEnergy

def make_rdc_dictionary(energy_vs_rdc):
    energy_vs_rdc_dict = {"energy_mev":list(energy_vs_rdc[:,0]), "rdc":energy_vs_rdc[:,1]}

def reference_particle_fit(qual_data, reference_energy, fit_type=None, fit_parameters=None, zero_fit=None):
    """

    Args:
        qual_data:
        reference_energy:
        fit_type:
        fit_parameters:
        zero_fit:

    Returns:

    """
    particle_energies = np.unique(qual_data[:, 0])
    index_of_energy_to_normalize_rdc = list(particle_energies).index(reference_energy)#np.argwhere(particle_energies == reference_energy)
    qual_data_grouped_by_particle_energy = group_qual_data_by_particle_energy(qual_data)

    if zero_fit == 'y':
        qual_data = qual_data_grouped_by_particle_energy[index_of_energy_to_normalize_rdc][:, [1, 2]]
        qual_data = np.vstack(([0,1],qual_data))

    else:
        qual_data = qual_data_grouped_by_particle_energy[index_of_energy_to_normalize_rdc][:, [1, 2]]

    if fit_type.lower() == 'single':
        reference_electron_fit_object = d_eq.singleDegradationEquation(qual_data)
        reference_electron_fit_object.fit(fit_parameters)

    elif fit_type.lower() == 'double':
        reference_electron_fit_object = d_eq.doubleDegradationEquation(qual_data)
        reference_electron_fit_object.fit(fit_parameters)

    else:
        reference_electron_fit_object = d_eq.singleDegradationEquation(qual_data)
        reference_electron_fit_object.fit(fit_parameters)

    return reference_electron_fit_object

def qual_data_from_file(rdc_txt_file):
    """
    Gets rdc curve from txt file from qual report and interpolates it
    Args:
        rdc_txt_file:

    Returns:
        interpolated rdc curve
    """
    qual_data_rdc_file = rdc_txt_file
    qual_data_rdc_from_file = np.loadtxt(qual_data_rdc_file, delimiter='\t', skiprows=1)
    intp = sp.interpolate.interp1d(qual_data_rdc_from_file[:, 0], qual_data_rdc_from_file[:, 1],
                                   fill_value='extrapolate')
    rdc_x_range = np.linspace(np.min(qual_data_rdc_from_file[:,0]), np.max(qual_data_rdc_from_file[:,0]),3000)
    new_rdc = np.zeros((len(rdc_x_range),2))
    new_rdc[:,0] = np.copy(rdc_x_range)
    new_rdc[:, 1] = intp(rdc_x_range)
    return new_rdc

# def plot_rad_data(qual_data, fit_type='Single', fit_parameters=None, save_fig='n', fig_name='rad fits'):
#     """
#     Checks to make sure fits of radiation data are good
#     Args:
#         qual_data:
#         param:
#         fit_type:
#         fit_parameters:
#         save_fig:
#         fig_name:
#
#     Returns:
#
#     """
#     print('rad fits')
#     plt.figure(fig_name)
#     particle_energy = np.unique(getattr(qual_data, param)[:,0])
#     groups = qual_data.group_by_energy(param)
#     colors = ['blue', 'green', 'orange', 'purple', 'pink', 'red', 'cyan', 'black']
#     min_fluence = []
#     max_fluence = []
#     for i, group in enumerate(groups):
#         min_fluence.append(group[0, 1])
#         max_fluence.append(group[-1, 1])
#         # group = np.vstack((np.array([[group[0,0],0,1]]), group))
#         plt.plot(group[:,1], group[:,2], 'o', color=colors[i], label=str(particle_energy[i])+' MeV')
#         # print(group[:, [1,2]])
#         if fit_type == 'Single':
#             rdc_model = d_eq.singleDegradationEquation(group[:, [1,2]])
#             rdc_model.fit(parameters=None)
#         elif fit_type =='Double':
#             rdc_model = d_eq.doubleDegradationEquation(group[:, [1, 2]])
#             rdc_model.fit(parameters=fit_parameters)
#         print(' {} MeV r2 rdc fits: {}, coef: {}'.format(group[0,0], rdc_model.r2, rdc_model.coefficients))
#
#         x = np.logspace(np.log10(group[0,1]), np.log10(group[-1,1]), 100)
#         # x = np.logspace(np.log10(group[0, 1])-0.2, np.log10(group[-1, 1]), 100)
#         plt.plot(x, rdc_model.getRemainingFactor(x), '--', color=colors[i])
#
#     # plt.xscale('symlog', linthreshx = 1e9)
#     plt.xscale('log')
#     plt.legend()
#     plt.xlabel(r'Fluence (p/cm$^2$)')
#     plt.ylabel(r'param')
