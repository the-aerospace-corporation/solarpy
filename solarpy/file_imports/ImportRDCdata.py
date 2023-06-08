import os
import pandas as pd
import numpy as np

class ImportRDCdata:
    def __init__(self, rdc_data_xlsx):
        self.filename = os.path.basename(rdc_data_xlsx)
        self.qualDataFrame = pd.ExcelFile(rdc_data_xlsx)
        self.particle_types = self.qualDataFrame.sheet_names
        self.electron_rdc_interpolated = []
        self.proton_rdc_interpolated = []
        self.proton_to_electron_conv = []
        self.shield_thickness_cm = []
        for particle_type in self.qualDataFrame.sheet_names:
            if particle_type in ['Electron RDCs Interpolated']:
                self.electron_rdc_interpolated = parseRadiationQualificationData(self.qualDataFrame.parse(particle_type))
            elif particle_type in ['Proton RDCs Interpolated']:
                self.proton_rdc_interpolated = parseRadiationQualificationData(self.qualDataFrame.parse(particle_type))
            elif  particle_type in ['Cell RDC Properties']:
                self.proton_to_electron_conv = self.qualDataFrame.parse(particle_type).values.astype('float')[0]
                self.shield_thickness_cm = self.qualDataFrame.parse(particle_type).values.astype('float')[1]



class parseRadiationQualificationData:
    def __init__(self, qualData):
        self.qualData = qualData
        self.data = qualData.values.astype('float')
        self.rdcs_available = list(qualData.columns[1:])
        self.particleEnergy = self.data[:,0]
        self.Voc = []
        self.Jsc = []
        self.Vmax = []
        self.Isc = []
        self.Imax = []
        self.Pmax = []
        self.FillFactor = []
        self.Efficiency = []
        self.remainingFactorsColumnIndices = []

        for i, column in enumerate(self.qualData.columns):
            if ('Voc' in column):
                if ~np.isnan(self.data[0,i]):
                    self.Voc = self.data[:,[0,i]]
            elif ('Jsc' in column):
                if ~np.isnan(self.data[0,i]):
                    self.Jsc = self.data[:,[0,i]]
            elif ('Vmax' in column):
                if ~np.isnan(self.data[0,i]):
                    self.Vmax = self.data[:,[0,i]]
            elif ('Isc' in column):
                if ~np.isnan(self.data[0,i]):
                    self.Isc = self.data[:,[0,i]]
            elif ('Imax' in column):
                if ~np.isnan(self.data[0,i]):
                    self.Imax = self.data[:,[0,i]]
            elif ('Pmax' in column):
                if ~np.isnan(self.data[0,i]):
                    self.Pmax = self.data[:,[0,i]]
                    # self.Efficiency = self.Pmax
            elif (column in ['FillFactor', 'FF', 'Fill Factor']):
                if self.data[0,i] and ~np.isnan(self.data[0,i]):
                    self.FillFactor = self.data[:,[0,i]]
                    self.rdcs_available[i-1] = u'FillFactor'  # minus 2 is because remaining Factor titles start after particle energy and fluence columns
            elif (column in ['Efficiency', 'eff', 'Eff', 'EFF']):
                if self.data[0,i] and ~np.isnan(self.data[0,i]):
                    self.Efficiency = self.data[:,[0,i]]
                    # self.Pmax = self.Efficiency
                    self.rdcs_available[i - 1] = u'Efficiency'

    def _groupByEnergy(self):
        indicesOfDuplicates = []
        energies = np.unique(self.data[:,0])
        for energy in energies:
            indicesOfDuplicates.append(list(np.where(self.data[:,0] == energy)[0]))
        groupedByEnergy = []
        for index in indicesOfDuplicates:
            groupedByEnergy.append(self.data[min(index):max(index) + 1])
        return groupedByEnergy
