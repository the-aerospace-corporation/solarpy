import numpy as np
from solarpy import degradation_equations as d_eq
from solarpy.eqflux import electron_rdc
from solarpy import data
from solarpy.eqflux.electron_rdc import electron_rdc_aero
from solarpy.eqflux.proton_rdc import proton_rdc_aero
from solarpy.examples.EQFLUX_electrons import EQFLUX_electrons
from solarpy.examples.EQFLUX_protons import EQFLUX_protons
import solarpy.eqflux.relativeMeV_fluence as rel_MeV
from solarpy.analyze_goodness_of_fit import analyze_goodness_of_fit


class eqflux_aero(object):
    """
    Class that derives all EQFLUX parameters for a solar cell including, the RDC of electrons, protons as well as the remaining factor.  If giving the environment spectrum the remaining factor of the solar cell can be calculated.
    
    Args:
        electron_qual_data (ndarray): [description]
        proton_qual_data (ndarray): [description]
        electron_differential_spectrum (ndarray): [description]
        proton_differential_spectrum (ndarray): [description]
        electron_reference_energy (float): Defaults to None. [description]
        proton_reference_energy (float): Defaults to None. [description]
        fit_type (str): Defaults to None. [description]
        critical_remaining_factor (float): Defaults to None. [description]
        density_of_shield (float): Defaults to None. [description]
        shield_thickness_cm (float): Defaults to None. [description]
        proton_shield_range_table (ndarray): Defaults to None. [description]
        electron_shield_range_table (ndarray): Defaults to None. [description]
        solar_cell_range_table (ndarray): Defaults to None. [description]
        solar_cell_displacement_table (ndarray): Defaults to None. [description]
        solar_cell_thickness_cm (ndarray): Defaults to None. [description]
        electron_energy_vs_rdc (ndarray): Defaults to None. [description]
        proton_energy_vs_rdc (ndarrayl): Defaults to None. [description]
        proton_to_electron_conversion_factor (float): Defaults to None. [description]
    """
    def __init__(self, electron_qual_data, proton_qual_data, electron_differential_spectrum, proton_differential_spectrum, electron_reference_energy=None, proton_reference_energy=None, critical_remaining_factor=None, density_of_shield=None, shield_thickness_cm=None, proton_shield_range_table=None, electron_shield_range_table=None, solar_cell_range_table=None,
                 solar_cell_displacement_table=None, solar_cell_thickness_cm=None, electron_energy_vs_rdc=None, proton_energy_vs_rdc=None, proton_to_electron_conversion_factor=None, p_fit_type=None, p_fit_parameters=None, e_fit_type=None, e_fit_parameters=None):


        self.electron_qual_data = electron_qual_data
        self.proton_qual_data = proton_qual_data
        self.electron_differential_spectrum = electron_differential_spectrum
        self.proton_differential_spectrum = proton_differential_spectrum

        if electron_reference_energy is None:
            electron_reference_energy = 1.0

        if proton_reference_energy is None:
            proton_reference_energy = 10

        if critical_remaining_factor is None:
            critical_remaining_factor = 0.8 # GaAs Radiation Handbook value

        if density_of_shield is None:
            density_of_shield = 2.554 # default to CMG glass density

        if shield_thickness_cm is None:
            shield_thickness_cm = 0 # default to no shield

        if electron_shield_range_table is None:
            electron_shield_range_table = data.CMG_electron_ranges

        if proton_shield_range_table is None:
            proton_shield_range_table = data.CMG_proton_ranges

        if solar_cell_range_table is None:
            solar_cell_range_table = data.GaAs_proton_ranges

        if solar_cell_displacement_table is None:
            solar_cell_displacement_table = data.GaAs_total_displacements

        if solar_cell_thickness_cm is None:
            solar_cell_thickness_cm = 5e-3

        self.electron_reference_energy = electron_reference_energy
        self.proton_reference_energy = proton_reference_energy
        self.e_fit_type = e_fit_type
        self.p_fit_type = p_fit_type
        self.e_fit_parameters = e_fit_parameters
        self.p_fit_parameters = p_fit_parameters
        self.critical_remaining_factor = critical_remaining_factor
        self.density_of_shield = density_of_shield
        self.shield_thickness_cm = shield_thickness_cm
        self.electron_shield_range_table = electron_shield_range_table
        self.proton_shield_range_table = proton_shield_range_table
        self.solar_cell_range_table = solar_cell_range_table
        self.solar_cell_displacement_table = solar_cell_displacement_table
        self.solar_cell_thickness_cm = solar_cell_thickness_cm

        self.reference_electron_energy_object = self.reference_particle_fit(self.electron_qual_data, self.electron_reference_energy, self.e_fit_type, self.e_fit_parameters)
        self.reference_proton_energy_object = self.reference_particle_fit(self.proton_qual_data, self.proton_reference_energy, self.p_fit_type, self.p_fit_parameters)

        self.electron_rdc_object = electron_rdc_aero(self.electron_qual_data, self.critical_remaining_factor, self.electron_reference_energy, self.density_of_shield, self.shield_thickness_cm, self.electron_shield_range_table, self.e_fit_type, self.e_fit_parameters)
        self.proton_rdc_object = proton_rdc_aero(self.proton_qual_data, self.critical_remaining_factor, self.proton_reference_energy, self.density_of_shield, self.shield_thickness_cm, self.proton_shield_range_table, self.solar_cell_range_table, self.solar_cell_displacement_table, self.solar_cell_thickness_cm, self.p_fit_type, self.p_fit_parameters)

        if electron_energy_vs_rdc is None:
            electron_energy_vs_rdc = self.electron_rdc_object.omnidirectional_shielded_rdc

        if proton_energy_vs_rdc is None:
            proton_energy_vs_rdc = self.proton_rdc_object.omnidirectional_shielded_rdc

        if proton_to_electron_conversion_factor is None:
            proton_to_electron_conversion_factor = self.get_proton_to_electron_factor()

        self.proton_to_electron_conversion_factor = proton_to_electron_conversion_factor
        self.electron_energy_vs_rdc = electron_energy_vs_rdc
        self.proton_energy_vs_rdc = proton_energy_vs_rdc

        # self.relative_electrons = EQFLUX_electrons(self.electron_energy_vs_rdc, self.electron_differential_spectrum)
        # self.relative_protons = EQFLUX_protons(self.proton_energy_vs_rdc, self.proton_differential_spectrum)

        self.one_MeV_electrons = []
        self.ten_MeV_protons = []
        self.one_MeV_electrons_from_10_MeV_proton = []
        self.total_one_MeV_electrons = []

        self.update()

    def update(self):
        """
        Updates member variables
        """
        intp_electron_differential_spectrum = rel_MeV.interpolate_electron_spectrum(self.electron_differential_spectrum)
        intp_proton_differential_spectrum = rel_MeV.interpolate_proton_spectrum(self.proton_differential_spectrum)
        # self.relative_electrons = EQFLUX_electrons(self.electron_energy_vs_rdc, self.electron_differential_spectrum)
        # self.relative_protons = EQFLUX_protons(self.proton_energy_vs_rdc, self.proton_differential_spectrum)
        self.one_MeV_electrons = rel_MeV.get_relative_MeV_fluence_trapezoidal_integration(self.electron_energy_vs_rdc, intp_electron_differential_spectrum)#self.relative_electrons.rel_MeV_Fluence()
        self.ten_MeV_protons = rel_MeV.get_relative_MeV_fluence_trapezoidal_integration(self.proton_rdc_object, intp_proton_differential_spectrum)
        self.one_MeV_electrons_from_10_MeV_proton = self.ten_MeV_protons * self.proton_to_electron_conversion_factor
        self.total_one_MeV_electrons = self.one_MeV_electrons_from_10_MeV_proton + self.one_MeV_electrons
        self.remaining_factor = self.reference_electron_energy_object.getRemainingFactor(self.total_one_MeV_electrons)

    def update_rdc(self):
        self.electron_rdc_object = electron_rdc_aero(self.electron_qual_data, self.critical_remaining_factor, self.electron_reference_energy, self.density_of_shield, self.shield_thickness_cm, self.electron_shield_range_table, self.e_fit_type, self.e_fit_parameters)
        self.proton_rdc_object = proton_rdc_aero(self.proton_qual_data, self.critical_remaining_factor, self.proton_reference_energy, self.density_of_shield, self.shield_thickness_cm, self.proton_shield_range_table, self.solar_cell_range_table, self.solar_cell_displacement_table, self.solar_cell_thickness_cm, self.p_fit_type, self.p_fit_parameters)
        self.electron_energy_vs_rdc = self.electron_rdc_object.omnidirectional_shielded_rdc
        self.proton_energy_vs_rdc = self.proton_rdc_object.omnidirectional_shielded_rdc

    def get_proton_to_electron_factor(self):
        proton_remaining_factor = self.reference_proton_energy_object.getFluence(self.critical_remaining_factor)
        electron_remaining_factor = self.reference_electron_energy_object.getFluence(self.critical_remaining_factor)
        proton_to_electron_factor = electron_remaining_factor/proton_remaining_factor
        return proton_to_electron_factor

    def reference_particle_fit(self, qual_data, reference_energy, fit_type=None, fit_parameters=None):
        return electron_rdc.reference_particle_fit(qual_data, reference_energy,fit_type, fit_parameters)

    def adjust_electron_rdc_fit_type(self, particle=None, fit_type=None, fit_parameters=None):
        if particle == 'p':
            rdc_object = self.proton_rdc_object
        elif particle == 'e':
            rdc_object = self.electron_rdc_object

        if fit_type == 'Single':
            rdc_object.get_energy_vs_rdc(fit_type=fit_type, fit_parameters=fit_parameters)
            rdc_object.get_energy_vs_rdc_extrapolated()
            rdc_object.get_omnidirectional_shielded_rdc()
            if particle == 'e':
                self.electron_energy_vs_rdc = rdc_object.omnidirectional_shielded_rdc
            elif particle == 'p':
                self.proton_energy_vs_rdc = rdc_object.omnidirectional_shielded_rdc

        elif fit_type == 'Double':
            rdc_object.get_energy_vs_rdc(fit_type=fit_type, fit_parameters=fit_parameters)
            rdc_object.get_energy_vs_rdc_extrapolated()
            rdc_object.get_omnidirectional_shielded_rdc()
            if particle == 'e':
                self.electron_energy_vs_rdc = rdc_object.omnidirectional_shielded_rdc
            elif particle == 'p':
                self.proton_energy_vs_rdc = rdc_object.omnidirectional_shielded_rdc

    def get_r2(self, particle,fit_type=None, fit_parameters=None):
        if particle == 'e':
            qual_data = self.electron_qual_data
        elif particle == 'p':
            qual_data = self.proton_qual_data

        qualDataGroupedByParticleEnergy = electron_rdc.group_qual_data_by_particle_energy(qual_data)

        for i, qual_data_by_energy in enumerate(qualDataGroupedByParticleEnergy):
            fluence_vs_remaining_factor = qual_data_by_energy[:, [1, 2]]
            if fit_type=='Single':
                model = d_eq.singleDegradationEquation(fluence_vs_remaining_factor)
                model.fit(fit_parameters)

            elif fit_type=='Double':
                model = d_eq.doubleDegradationEquation(fluence_vs_remaining_factor)
                model.fit(fit_parameters)

            else:
                model = d_eq.singleDegradationEquation(fluence_vs_remaining_factor)
                model.fit(fit_parameters)
            print('{} MeV r2: {}, coef: {}'.format(qual_data_by_energy[0,0], model.r2, model.coefficients))
