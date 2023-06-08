import os
import pandas as pd

class ImportRangeTable:
    def __init__(self, qual_data_xlsx):
        self.filename = os.path.basename(qual_data_xlsx)
        self.qualDataFrame = pd.ExcelFile(qual_data_xlsx)
        self.particle_types = self.qualDataFrame.sheet_names
        self.electron_range_table = []
        self.proton_range_table = []

        for particle_type in self.qualDataFrame.sheet_names:
            if particle_type in ['Electrons', 'electrons', 'Electron', 'electron']:
                self.electron_range_table = self.qualDataFrame.parse(particle_type)
            elif particle_type in ['Protons', 'protons', 'Proton', 'proton']:
                self.proton_range_table = self.qualDataFrame.parse(particle_type)
