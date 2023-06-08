from solarpy.file_imports.ImportRadiationQualData import ImportRadiationQualData
from solarpy.CLI.ElectronRDC import ElectronRDC

class SolarCellElectronRDC(ImportRadiationQualData):
    def __init__(self, qual_data_xlsx, critical_remaining_factor=None, energy_to_normalize_rdc=None, fit_type=None, fit_parameters=None, density_of_shield=None,
                 shield_thickness_cm=None, shield_range_table=None):
        ImportRadiationQualData.__init__(self, qual_data_xlsx)
        self.remainingFactorsAvailable = self.electronQualData.remainingFactorsAvailable
        self.Voc = []
        self.Jsc = []
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
        self.electron_shield_range_table = shield_range_table

        for i, data_rf in enumerate(self.remainingFactorsAvailable):
            if data_rf == 'Voc':
                if getattr(self.electronQualData, data_rf)[0,0]:
                    self.Voc = self.electron_rdc(getattr(self.electronQualData, data_rf))
            elif data_rf == 'Isc':
                if getattr(self.electronQualData, data_rf)[0,0]:
                    self.Isc = self.electron_rdc(getattr(self.electronQualData, data_rf))
            elif data_rf == 'Vmax':
                if getattr(self.electronQualData, data_rf)[0,0]:
                    self.Vmax = self.electron_rdc(getattr(self.electronQualData, data_rf))
            elif data_rf == 'Imax':
                if getattr(self.electronQualData, data_rf)[0,0]:
                    self.Imax = self.electron_rdc(getattr(self.electronQualData, data_rf))
            elif data_rf == 'FillFactor':
                if getattr(self.electronQualData, data_rf)[0,0]:
                    self.FillFactor = self.electron_rdc(getattr(self.electronQualData, data_rf))
            elif data_rf == 'Pmax':
                if getattr(self.electronQualData, data_rf)[0,0]:
                    self.Pmax = self.electron_rdc(getattr(self.electronQualData, data_rf))
            elif data_rf == 'Efficiency':
                if getattr(self.electronQualData, data_rf)[0,0]:
                    self.Efficiency = self.electron_rdc(getattr(self.electronQualData, data_rf))
            else:
                print('no parameters found')

    def electron_rdc(self, qual_data):
        rdc_object = ElectronRDC(qual_data,
                                 critical_remaining_factor=self.critical_remaining_factor,
                                 energy_to_normalize_rdc=self.energy_to_normalize_rdc,
                                 fit_type=self.fit_type,
                                 fit_parameters=self.fit_parameters,
                                 density_of_shield=self.density_of_shield,
                                 shield_thickness_cm=self.shield_thickness_cm,
                                 shield_range_table=self.electron_shield_range_table
                                 )
        return rdc_object

