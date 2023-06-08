from solarpy.file_imports.ImportRadiationQualData import ImportRadiationQualData
from solarpy.CLI.ProtonRDC import ProtonRDC
class SolarCellProtonRDC(ImportRadiationQualData):
    def __init__(self, qual_data_xlsx, critical_remaining_factor=None, energy_to_normalize_rdc=None, fit_type=None, fit_parameters=None, density_of_shield=None,
                 shield_thickness_cm=None, shield_range_table=None, solar_cell_range_table=None,
                 solar_cell_displacement_table=None, solar_cell_thickness_cm=None):
        ImportRadiationQualData.__init__(self, qual_data_xlsx)
        self.remainingFactorsAvailable = self.protonQualData.remainingFactorsAvailable
        self.Voc = []
        self.Vmax = []
        self.Isc = []
        self.Imax = []
        self.Pmax = []
        self.FillFactor = []
        self.Efficiency = []

        self.critical_remaining_factor = critical_remaining_factor
        self.energy_to_normalize_rdc = energy_to_normalize_rdc
        self.fit_type = fit_type
        self.fit_parameters = fit_parameters
        self.density_of_shield = density_of_shield
        self.shield_thickness_cm = shield_thickness_cm
        self.proton_shield_range_table = shield_range_table
        self.solar_cell_range_table = solar_cell_range_table
        self.solar_cell_displacement_table = solar_cell_displacement_table
        self.solar_cell_thickness_cm = solar_cell_thickness_cm
        for i, data_rf in enumerate(self.remainingFactorsAvailable):
            if data_rf == 'Voc':
                if getattr(self.protonQualData, data_rf)[0,0]:
                    self.Voc = self.proton_rdc(getattr(self.protonQualData, data_rf))
            elif data_rf == 'Isc':
                if getattr(self.protonQualData, data_rf)[0,0]:
                    self.Isc = self.proton_rdc(getattr(self.protonQualData, data_rf))
            elif data_rf == 'Vmax':
                if getattr(self.protonQualData, data_rf)[0,0]:
                    self.Vmax = self.proton_rdc(getattr(self.protonQualData, data_rf))
            elif data_rf == 'Imax':
                if getattr(self.protonQualData, data_rf)[0,0]:
                    self.Imax = self.proton_rdc(getattr(self.protonQualData, data_rf))
            elif data_rf == 'FillFactor':
                if getattr(self.protonQualData, data_rf)[0,0]:
                    self.FillFactor = self.proton_rdc(getattr(self.protonQualData, data_rf))
            elif data_rf == 'Pmax':
                if getattr(self.protonQualData, data_rf)[0,0]:
                    self.Pmax = self.proton_rdc(getattr(self.protonQualData, data_rf))
            elif data_rf == 'Efficiency':
                if getattr(self.protonQualData, data_rf)[0,0]:
                    self.Efficiency = self.proton_rdc(getattr(self.protonQualData, data_rf))
            else:
                print('no parameters found')

    def proton_rdc(self, qual_data):
        rdc_object = ProtonRDC(qual_data,
                               critical_remaining_factor=self.critical_remaining_factor,
                               energy_to_normalize_rdc=self.energy_to_normalize_rdc,
                               fit_type=self.fit_type,
                               fit_parameters=self.fit_parameters,
                               density_of_shield=self.density_of_shield,
                               shield_thickness_cm=self.shield_thickness_cm,
                               shield_range_table=self.proton_shield_range_table,
                               solar_cell_range_table= self.solar_cell_range_table,
                               solar_cell_displacement_table=self.solar_cell_displacement_table,
                               solar_cell_thickness_cm=self.solar_cell_thickness_cm
                               )
        return rdc_object
