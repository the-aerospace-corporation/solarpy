import pandas as pd

# class OutputRDC_


class OutputCellTechRDCdata(object):
    def __init__(self, electron_rdc=None, proton_rdc=None, type_of_rdc=None, params=None ):
        """
        Used for CLI
        Args:
            electron_rdc: object from SolarCellElectronRDC
            proton_rdc: object from SolarCellProtonRDC
            type_of_rdc:
            params:
        """

        self.rfs = ['Voc', 'Isc', 'Vmax', 'Imax', 'FillFactor', 'Pmax', 'Efficiency']
        if (params is None) or (params == 'All'):
            self.params = self.rfs
        else:
            self.params = params
        self.electron_rdc_dataframe = pd.DataFrame(columns=['Energy (MeV)'] + self.rfs)
        self.proton_rdc_dataframe = pd.DataFrame(columns=['Energy (MeV)'] + self.rfs)
        self.type_of_rdc = type_of_rdc

        if electron_rdc:
            self.fill_dataframe(electron_rdc, self.electron_rdc_dataframe)

        if electron_rdc:
            self.fill_dataframe(proton_rdc, self.proton_rdc_dataframe)


    def fill_dataframe(self, rdc, dataframe):
        for param in self.params:
            if param in rdc.remainingFactorsAvailable:
                if self.type_of_rdc in ['unidirectional', 'u']:
                    rdc_curve = getattr(rdc, param).get_energy_vs_rdc()

                elif self.type_of_rdc in ['unidirectional extrapolated', 'ui']:
                    rdc_curve = getattr(rdc, param).get_energy_vs_rdc_extrapolated()

                elif self.type_of_rdc in ['omnidirectional shielded', 'os']:
                    rdc_curve = getattr(rdc, param).get_omnidirectional_shielded_rdc()

                else:
                    rdc_curve = getattr(rdc, param).get_energy_vs_rdc_extrapolated()

                dataframe['Energy (MeV)'] = rdc_curve[:, 0]
                dataframe[param] = rdc_curve[:, 1]

    def output(self, filename):
        filename = filename+'.xlsx'
        with pd.ExcelWriter(filename) as writer:
            self.electron_rdc_dataframe.to_excel(writer, sheet_name='Electron RDCs', index=False)
            self.proton_rdc_dataframe.to_excel(writer, sheet_name='Proton RDCs', index=False)


class OutputCellTechQualData(object):
    def __init__(self,electron_qual_data=None, proton_qual_data=None):
        self.rfs = ['Voc', 'Isc', 'Vmax', 'Imax', 'FillFactor', 'Pmax', 'Efficiency']
        self.electron_qual_dataframe = pd.DataFrame(columns=['Energy','Fluence(e/cm2)'] + self.rfs)
        self.proton_qual_dataframe = pd.DataFrame(columns=['Energy','Fluence(p/cm2)'] + self.rfs)
        self.electron_rdc_dataframe = pd.DataFrame(columns=['Energy', 'RDC'])
        self.proton_rdc_dataframe = pd.DataFrame(columns=['Energy', 'RDC'])

        if electron_qual_data:
            self.electron_qual_dataframe['Energy (MeV)'] = electron_qual_data.particleEnergy
            self.electron_qual_dataframe['Fluence(e/cm2)'] = electron_qual_data.fluence
            self.fill_dataframe(electron_qual_data, self.electron_qual_dataframe)
            # self.electron_qual_dataframe.fillna(0, inplace=True)

        if proton_qual_data:
            self.proton_qual_dataframe['Energy (MeV)'] = proton_qual_data.particleEnergy
            self.proton_qual_dataframe['Fluence(p/cm2)'] = proton_qual_data.fluence
            self.fill_dataframe(proton_qual_data, self.proton_qual_dataframe)
            # self.proton_qual_dataframe.fillna(0, inplace=True)

    def fill_dataframe(self, qual_data, dataframe):
        for rf in self.rfs:
            if rf in qual_data.remainingFactorsAvailable:
                dataframe[rf] = getattr(qual_data, rf)[:,-1]

    def output(self, filename):
        filename = filename+'.xlsx'
        with pd.ExcelWriter(filename) as writer:
            self.electron_qual_dataframe.to_excel(writer, sheet_name='Electron Rad Data', index=False)
            self.proton_qual_dataframe.to_excel(writer, sheet_name='Proton Rad Data', index=False)
