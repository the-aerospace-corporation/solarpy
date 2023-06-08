import pandas as pd
import json


class OutputCellTechRDCdata(object):
    def __init__(self, electron_rdc=None, proton_rdc=None, type_of_rdc=None, params=None):
        """
        Used for CLI
        Args:
            electron_rdc: object from SolarCellElectronRDC
            proton_rdc: object from SolarCellProtonRDC
            type_of_rdc:
            params:
        """
        self.electron_rdc_object = electron_rdc
        self.proton_rdc_object = proton_rdc
        self.rfs = ['Voc', 'Isc', 'Vmax', 'Imax', 'FillFactor', 'Pmax', 'Efficiency']
        self.rfs_json = ['v_oc', 'i_sc', 'v_max', 'i_max', 'fill_factor', 'p_max', 'efficiency']

        if (params is None) or (params == 'All'):
            self.params = self.rfs
        else:
            self.params = params
        self.type_of_rdc = type_of_rdc

    def fill_dataframe(self, rdc, type_of_rdc):
        rdc_dataframe = pd.DataFrame(columns=['Energy (MeV)'] + self.rfs)
        rdc_dataframe_json = pd.DataFrame(columns=['energy_mev)'] + self.rfs_json)

        for param in self.params:
            if param in rdc.remainingFactorsAvailable:
                if getattr(rdc, param):

                    if type_of_rdc in ['unidirectional', 'u']:
                        rdc_curve = getattr(rdc, param).get_energy_vs_rdc()

                    elif type_of_rdc in ['unidirectional extrapolated', 'ui']:
                        rdc_curve = getattr(rdc, param).get_energy_vs_rdc_extrapolated()

                    elif type_of_rdc in ['omnidirectional shielded', 'os']:
                        rdc_curve = getattr(rdc, param).get_omnidirectional_shielded_rdc()

                    elif type_of_rdc in ['convert', 'c']:
                        rdc_curve = rdc
                    else:
                        rdc_curve = getattr(rdc, param).get_energy_vs_rdc_extrapolated()

                    rdc_dataframe['Energy (MeV)'] = rdc_curve[:, 0]
                    rdc_dataframe[param] = rdc_curve[:, 1]
                    # rdc_dataframe[]

                else:
                    # rdc_dataframe['Energy (MeV)'] = []
                    rdc_dataframe[param] = pd.Series([])

        return rdc_dataframe

    def output(self, filename):

        if self.electron_rdc_object:
            electron_rdc_dataframe = self.fill_dataframe(self.electron_rdc_object, self.type_of_rdc)

        if self.proton_rdc_object:
            proton_rdc_dataframe = self.fill_dataframe(self.proton_rdc_object, self.type_of_rdc)

        filename = filename + '.xlsx'
        with pd.ExcelWriter(filename) as writer:
            electron_rdc_dataframe.to_excel(writer, sheet_name='Electron RDCs', index=False)
            proton_rdc_dataframe.to_excel(writer, sheet_name='Proton RDCs', index=False)

    def output_json(self, filename):
        if self.electron_rdc_object:
            electron_rdc_dataframe = self.fill_dataframe(self.electron_rdc_object, self.type_of_rdc).fillna(0, inplace=True)

        if self.proton_rdc_object:
            proton_rdc_dataframe = self.fill_dataframe(self.proton_rdc_object, self.type_of_rdc).fillna(0, inplace=True)

        filename_electron_rdc = filename + '_electron_rdc.json'
        electron_rdc_dataframe.to_json(filename_electron_rdc, orient='columns')

        filename_proton_rdc = filename + '_proton_rdc.json'
        proton_rdc_dataframe.to_json(filename_proton_rdc, orient='columns')

    def output_all(self, filename):
        if self.electron_rdc_object:
            electron_u_rdc_dataframe = self.fill_dataframe(self.electron_rdc_object, 'u')
            electron_ui_rdc_dataframe = self.fill_dataframe(self.electron_rdc_object, 'ui')
            electron_os_rdc_dataframe = self.fill_dataframe(self.electron_rdc_object, 'os')

        if self.proton_rdc_object:
            proton_u_rdc_dataframe = self.fill_dataframe(self.proton_rdc_object, 'u')
            proton_ui_rdc_dataframe = self.fill_dataframe(self.proton_rdc_object, 'ui')
            proton_os_rdc_dataframe = self.fill_dataframe(self.proton_rdc_object, 'os')

        filename_xlsx = filename + '.xlsx'
        with pd.ExcelWriter(filename_xlsx) as writer:
            electron_u_rdc_dataframe.to_excel(writer, sheet_name='Electron RDCs Unidirectional', index=False)
            proton_u_rdc_dataframe.to_excel(writer, sheet_name='Proton RDCs Unidirectional', index=False)

            electron_ui_rdc_dataframe.to_excel(writer, sheet_name='Electron RDCs Interpolated', index=False)
            proton_ui_rdc_dataframe.to_excel(writer, sheet_name='Proton RDCs Interpolated', index=False)

            electron_os_rdc_dataframe.to_excel(writer, sheet_name='Electron RDCs Shielded', index=False)
            proton_os_rdc_dataframe.to_excel(writer, sheet_name='Proton RDCs Shielded', index=False)

            proton_to_electron_conversion = 0
            shield_thickness_in_cm = 0
            rdc_cell_properties = pd.DataFrame([proton_to_electron_conversion, shield_thickness_in_cm],
                                               index=['Proton to Electron Conversion', 'Shield thickness in cm'])
            rdc_cell_properties.to_excel(writer, sheet_name='Cell RDC Properties', index=True)

        ###need to fill NaNs because Json does not like it, but for excel its nice because it makes empty tables
        electron_u_rdc_dataframe.fillna(0, inplace=True)
        electron_ui_rdc_dataframe.fillna(0, inplace=True)
        electron_os_rdc_dataframe.fillna(0, inplace=True)
        proton_u_rdc_dataframe.fillna(0, inplace=True)
        proton_ui_rdc_dataframe.fillna(0, inplace=True)
        proton_os_rdc_dataframe.fillna(0, inplace=True)

        json_build = {'electron_rdcs_unidirectional': electron_u_rdc_dataframe.to_dict(),
                      'electron_rdcs_interpolated': electron_ui_rdc_dataframe.to_dict(),
                      'electron_rdcs_shielded': electron_os_rdc_dataframe.to_dict(),
                      'proton_rdcs_unidirectional': proton_u_rdc_dataframe.to_dict(),
                      'proton_rdcs_interpolated': proton_ui_rdc_dataframe.to_dict(),
                      'proton_rdcs_shielded': proton_os_rdc_dataframe.to_dict(),
                      'cell_rdc_properties': rdc_cell_properties.to_dict()}

        with open(filename + '.json', 'w') as f:
            json.dump(json_build, f)
