import pandas as pd

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
