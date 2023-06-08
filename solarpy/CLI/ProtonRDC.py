from solarpy.eqflux import electron_rdc
from solarpy.eqflux import proton_rdc
import numpy as np
from solarpy import degradation_equations as d_eq

class ProtonRDC(object):
    def __init__(self, qual_data, critical_remaining_factor=None, energy_to_normalize_rdc=None, fit_type=None, fit_parameters=None,
                 density_of_shield=None, shield_thickness_cm=None, shield_range_table=None, solar_cell_range_table=None,
                 solar_cell_displacement_table=None, solar_cell_thickness_cm=None):
        if critical_remaining_factor == None:
            self.critical_remaining_factor = 0.8

        else:
            self.critical_remaining_factor = critical_remaining_factor

        if energy_to_normalize_rdc == None:
            self.energy_to_normalize_rdc = 10

        else:
            self.energy_to_normalize_rdc = energy_to_normalize_rdc

        self.qual_data = qual_data

        self.fit_type = fit_type
        self.fit_parameter = fit_parameters

        self.shield_thickness_cm = shield_thickness_cm
        self.density_of_shield = density_of_shield
        self.mass_thickness_of_shield = density_of_shield * shield_thickness_cm
        self.solar_cell_range_table = solar_cell_range_table
        self.solar_cell_displacement_table = solar_cell_displacement_table
        self.solar_cell_thickness_cm = solar_cell_thickness_cm
        self.shield_range_table = shield_range_table
        self.energy_vs_rdc = []
        self.energy_vs_rdc_extrapolated = []
        self.omnidirectional_shielded_rdc = []

        self.get_energy_vs_rdc(fit_parameters=self.fit_parameter, fit_type=self.fit_type)
        # self.get_energy_vs_rdc_extrapolated()
        # self.get_omnidirectional_shielded_rdc()

    def get_energy_vs_rdc(self, fit_parameters=None, fit_type=None):
        if self.fit_type:
            fit_type = self.fit_type

        if fit_type == 'Single':
            energy_vs_rdc = electron_rdc.get_energy_vs_rdc(self.qual_data, self.critical_remaining_factor,
                                                                          self.energy_to_normalize_rdc, fit_parameters=fit_parameters)

        elif fit_type == 'Double':
            energy_vs_rdc = electron_rdc.get_energy_vs_rdc_double_degradation(self.qual_data,
                                                                                             self.critical_remaining_factor,
                                                                                             self.energy_to_normalize_rdc,
                                                                                             fit_parameters=fit_parameters)

        else:
            energy_vs_rdc = electron_rdc.get_energy_vs_rdc(self.qual_data, self.critical_remaining_factor,
                                                                          self.energy_to_normalize_rdc, fit_parameters=fit_parameters)
        self.energy_vs_rdc = energy_vs_rdc
        return self.energy_vs_rdc

    def get_energy_vs_rdc_extrapolated(self, minimum_particle_energy=1e-5, maximum_particle_energy=1e2, minimum_energy_rdc=None, maximum_rdc_energy = None,
                                       number_of_points=3000):
        if maximum_rdc_energy is None:
            maximum_rdc_energy = ([50,0.61],[100,0.56])
        if minimum_energy_rdc is None:
            minimum_energy_rdc = [0.02, 0.1]
        # self.energy_vs_rdc_extrapolated = rdc.extrapolate_RDC(self.energy_vs_rdc, minimum_particle_energy,
        #                                                       maximum_particle_energy, minimum_energy_rdc, maximum_rdc_energy, number_of_points)

        if not self.energy_vs_rdc.any():
            self.get_energy_vs_rdc()

        self.energy_vs_rdc_extrapolated = electron_rdc.extrapolate_RDC_loglog(self.energy_vs_rdc, minimum_particle_energy,
                                                                                             maximum_particle_energy, minimum_energy_rdc, maximum_rdc_energy, number_of_points)

        return self.energy_vs_rdc_extrapolated

    def get_omnidirectional_shielded_rdc(self):
        if not self.energy_vs_rdc_extrapolated.any():
            self.get_energy_vs_rdc_extrapolated()

        omnidirectional_shielded_rdc = proton_rdc.get_omnidirectional_and_shielded_proton_rdc(
            self.energy_vs_rdc_extrapolated, self.solar_cell_range_table, self.shield_range_table,
            self.solar_cell_displacement_table, self.solar_cell_thickness_cm, self.shield_thickness_cm)
        self.omnidirectional_shielded_rdc = omnidirectional_shielded_rdc
        return self.omnidirectional_shielded_rdc

    def get_fluence_at_reference(self, critical_remaining_factor):
        return self.reference_particle_fit_object().getFluence(critical_remaining_factor)

    def reference_particle_fit_object(self):
        particle_energies = np.unique(self.qual_data[:, 0])
        index_of_energy_to_normalize_rdc = list(particle_energies).index(self.energy_to_normalize_rdc) #np.argwhere(particle_energies == reference_energy)
        qual_data_grouped_by_particle_energy = electron_rdc.group_qual_data_by_particle_energy(self.qual_data)

        if self.fit_type == 'single':
            reference_electron_fit_object = d_eq.singleDegradationEquation(qual_data_grouped_by_particle_energy[index_of_energy_to_normalize_rdc][:, [1, 2]])

        elif self.fit_type == 'double':
            reference_electron_fit_object = d_eq.doubleDegradationEquation(qual_data_grouped_by_particle_energy[index_of_energy_to_normalize_rdc][:, [1, 2]])

        else:
            reference_electron_fit_object = d_eq.singleDegradationEquation(qual_data_grouped_by_particle_energy[index_of_energy_to_normalize_rdc][:, [1, 2]])

        return reference_electron_fit_object

    # def plot_rdc(self):
    #     ax = plt.loglog(self.energy_vs_rdc[:, 0], self.energy_vs_rdc[:, 1])
    #     return ax
    #
    # def plot_omnidirectional_shielded_rdc(self):
    #     ax = plt.loglog(self.omnidirectional_shielded_rdc[:, 0], self.omnidirectional_shielded_rdc[:, 1])
    #     return ax