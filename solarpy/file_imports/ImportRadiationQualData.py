import pandas as pd
import os
import numpy as np
from solarpy.file_imports.ParseRadiationQualificationData import ParseRadiationQualificationData

class ImportRadiationQualData:
    def __init__(self, qual_data_xlsx):
        self.filename = os.path.basename(qual_data_xlsx)
        self.qualDataFrame = pd.ExcelFile(qual_data_xlsx)
        self.particleTypes = self.qualDataFrame.sheet_names
        self.protonQualData = []
        self.electronQualData = []

        for particleType in self.qualDataFrame.sheet_names:
            if particleType in ['Electrons', 'electrons', 'Electron', 'electron', 'Electron Rad Data']:
                self.electronQualData = ParseRadiationQualificationData(self.qualDataFrame.parse(particleType))
            elif particleType in ['Protons', 'protons', 'Proton', 'proton', 'Proton Rad Data']:
                self.protonQualData = ParseRadiationQualificationData(self.qualDataFrame.parse(particleType))


    def _groupByEnergy(self):
        indicesOfDuplicates = []
        energies = np.unique(self.data[:,0])
        for energy in energies:
            indicesOfDuplicates.append(list(np.where(self.data[:,0] == energy)[0]))
        groupedByEnergy = []
        for index in indicesOfDuplicates:
            groupedByEnergy.append(self.data[min(index):max(index) + 1])
        return groupedByEnergy
