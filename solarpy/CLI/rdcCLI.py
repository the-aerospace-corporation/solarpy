import numpy as np
from solarpy import data
from solarpy.CLI.ElectronRDC import ElectronRDC
from solarpy.CLI.ProtonRDC import ProtonRDC
from solarpy.file_imports.ImportRadiationQualData import ImportRadiationQualData
from solarpy.eqflux import electron_rdc
from solarpy import degradation_equations as d_eq


#### Not there yet, have to seperate proton and electon????


class rdc_all_params_CLI(ImportRadiationQualData):
    def __init__(self, qual_data_xlsx, electron_reference_energy=None, proton_reference_energy=None, fit_type=None, critical_remaining_factor_electrons=None, critical_remaining_factor_protons=None, density_of_shield=None, shield_thickness_cm=None, proton_shield_range_table=None, electron_shield_range_table=None, solar_cell_range_table=None,
                 solar_cell_displacement_table=None, solar_cell_thickness_cm=None):
        ImportRadiationQualData.__init__(self, qual_data_xlsx)

        self.Voc = []
        self.Jsc = []
        self.Vmax = []
        self.Isc = []
        self.Imax = []
        self.Pmax = []
        self.FillFactor = []
        self.Efficiency = []
        self.remaining_factors = ['Voc', 'Isc', 'Vmax', 'Imax', 'FillFactor', 'Pmax', 'Efficiency']

        self.electron_reference_energy = electron_reference_energy
        self.proton_reference_energy = proton_reference_energy
        self.fit_type = fit_type
        self.critical_remaining_factor_electrons = critical_remaining_factor_electrons
        self.critical_remaining_factor_protons = critical_remaining_factor_protons
        self.density_of_shield = density_of_shield
        self.shield_thickness_cm = shield_thickness_cm
        self.proton_shield_range_table = proton_shield_range_table
        self.electron_shield_range_table = electron_shield_range_table
        self.solar_cell_range_table = solar_cell_range_table
        self.solar_cell_displacement_table = solar_cell_displacement_table
        self.solar_cell_thickness_cm = solar_cell_thickness_cm


        for i, rf in enumerate(self.remaining_factors):
            if rf == 'Voc':
                self.Voc = rdc_single_param_CLI(electron_qual_data=self.electronQualData.Voc,
                                                proton_qual_data=self.protonQualData.Voc)

            elif rf == 'Isc':
                self.Isc = rdc_single_param_CLI(electron_qual_data=self.electronQualData.Isc,
                                                proton_qual_data=self.protonQualData.Isc)

            elif rf == 'Vmax':
                self.Vmax = rdc_single_param_CLI(electron_qual_data=self.electronQualData.Vmax,
                                                proton_qual_data=self.protonQualData.Vmax)

            elif rf == 'Imax':
                self.Imax = rdc_single_param_CLI(electron_qual_data=self.electronQualData.Imax,
                                                proton_qual_data=self.protonQualData.Imax)

            elif rf == 'FillFactor':
                self.FillFactor = rdc_single_param_CLI(electron_qual_data=self.electronQualData.FillFactor,
                                                proton_qual_data=self.protonQualData.FillFactor)

            elif rf == 'Pmax':
                self.Pmax = rdc_single_param_CLI(electron_qual_data=self.electronQualData.Pmax,
                                                proton_qual_data=self.protonQualData.Pmax)

            elif rf == 'Efficiency':
                self.Efficiency = rdc_single_param_CLI(electron_qual_data=self.electronQualData.Efficiency,
                                                proton_qual_data=self.protonQualData.Efficiency)

            else:
                print('no parameters found')

    def get_rdc_single_param(self, qual_data_electrons, qual_data_protons):
        rdc_single_param = rdc_single_param_CLI(electron_qual_data=qual_data_electrons,
                                                proton_qual_data=qual_data_protons,
                                                electron_reference_energy=self.electron_reference_energy,
                                                proton_reference_energy=self.proton_reference_energy,
                                                fit_type=self.fit_type,
                                                critical_remaining_factor_electrons=self.critical_remaining_factor_electrons,
                                                critical_remaining_factor_protons=self.critical_remaining_factor_protons,
                                                density_of_shield=self.density_of_shield,
                                                shield_thickness_cm=self.shield_thickness_cm,
                                                proton_shield_range_table=self.proton_shield_range_table,
                                                electron_shield_range_table=self.electron_shield_range_table,
                                                solar_cell_range_table=self.solar_cell_range_table,
                                                solar_cell_displacement_table=self.solar_cell_displacement_table,
                                                solar_cell_thickness_cm=self.solar_cell_thickness_cm
                                                )
        return  rdc_single_param

class rdc_single_param_CLI(object):
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
    def __init__(self, electron_qual_data, proton_qual_data, electron_reference_energy=None, proton_reference_energy=None, fit_type=None, critical_remaining_factor_electrons=None, critical_remaining_factor_protons=None, density_of_shield=None, shield_thickness_cm=None, proton_shield_range_table=None, electron_shield_range_table=None, solar_cell_range_table=None,
                 solar_cell_displacement_table=None, solar_cell_thickness_cm=None):


        self.electron_qual_data = electron_qual_data
        self.proton_qual_data = proton_qual_data
        self.proton_to_electron_conversion_factor = []

        if electron_reference_energy is None:
            electron_reference_energy = 1.0

        if proton_reference_energy is None:
            proton_reference_energy = 10

        if fit_type is None:
            fit_type = 'single' # GaAs Radiation Handbook uses a polynomial...but everyone else seems to use this

        if critical_remaining_factor_electrons is None:
            critical_remaining_factor_electrons = 0.8 # GaAs Radiation Handbook value

        if critical_remaining_factor_protons is None:
            critical_remaining_factor_protons = 0.8 # GaAs Radiation Handbook value

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
        self.fit_type = fit_type
        self.critical_remaining_factor_electrons = critical_remaining_factor_electrons
        self.critical_remaining_factor_protons = critical_remaining_factor_protons
        self.density_of_shield = density_of_shield
        self.shield_thickness_cm = shield_thickness_cm
        self.electron_shield_range_table = electron_shield_range_table
        self.proton_shield_range_table = proton_shield_range_table
        self.solar_cell_range_table = solar_cell_range_table
        self.solar_cell_displacement_table = solar_cell_displacement_table
        self.solar_cell_thickness_cm = solar_cell_thickness_cm

        if electron_qual_data & proton_qual_data:
            self.reference_electron_energy_object = self.reference_particle_fit(self.electron_qual_data, self.electron_reference_energy, self.fit_type)
            self.reference_proton_energy_object = self.reference_particle_fit(self.proton_qual_data, self.proton_reference_energy, self.fit_type)

            self.electron_rdc = ElectronRDC(self.electron_qual_data, self.critical_remaining_factor_electrons, self.electron_reference_energy, self.density_of_shield, self.shield_thickness_cm, self.electron_shield_range_table, fit_type=self.fit_type)
            self.proton_rdc = ProtonRDC(self.proton_qual_data, self.critical_remaining_factor_protons, self.proton_reference_energy, self.density_of_shield, self.shield_thickness_cm, self.proton_shield_range_table, self.solar_cell_range_table, self.solar_cell_displacement_table, self.solar_cell_thickness_cm, fit_type=self.fit_type)

            self.proton_to_electron_conversion_factor = self.get_proton_to_electron_factor()

        elif electron_qual_data:
            self.electron_rdc = ElectronRDC(self.electron_qual_data, self.critical_remaining_factor_electrons, self.electron_reference_energy, self.density_of_shield, self.shield_thickness_cm, self.electron_shield_range_table)

        elif proton_qual_data:
            self.proton_rdc = ProtonRDC(self.proton_qual_data, self.critical_remaining_factor_protons, self.proton_reference_energy, self.density_of_shield, self.shield_thickness_cm, self.proton_shield_range_table, self.solar_cell_range_table, self.solar_cell_displacement_table, self.solar_cell_thickness_cm)

    def get_proton_to_electron_factor(self):
        proton_remaining_factor = self.reference_proton_energy_object.getFluence(self.critical_remaining_factor_protons)
        electron_remaining_factor = self.reference_electron_energy_object.getFluence(self.critical_remaining_factor_electrons)
        proton_to_electron_factor = electron_remaining_factor/proton_remaining_factor
        return proton_to_electron_factor

    def reference_particle_fit(self, qual_data, reference_energy, fit_type):
        particle_energies = np.unique(qual_data[:, 0])
        index_of_energy_to_normalize_rdc = list(particle_energies).index(reference_energy)#np.argwhere(particle_energies == reference_energy)
        qual_data_grouped_by_particle_energy = electron_rdc.group_qual_data_by_particle_energy(qual_data)

        if fit_type == 'single':
            reference_electron_fit_object = d_eq.singleDegradationEquation(qual_data_grouped_by_particle_energy[index_of_energy_to_normalize_rdc][:, [1, 2]])

        elif fit_type == 'double':
            reference_electron_fit_object = d_eq.doubleDegradationEquation(qual_data_grouped_by_particle_energy[index_of_energy_to_normalize_rdc][:, [1, 2]])

        else:
            reference_electron_fit_object = d_eq.singleDegradationEquation(qual_data_grouped_by_particle_energy[index_of_energy_to_normalize_rdc][:, [1, 2]])

        return reference_electron_fit_object
