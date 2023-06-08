import numpy as np

class ParseRadiationQualificationData:
    def __init__(self, qualData):
        self.qual_data = qualData
        self.data = qualData.values.astype('float')
        self.remaining_factors_available = list(qualData.columns[2:])
        self.particle_energy = self.data[:, 0]
        self.fluence = self.data[:,1]
        self.voc = []
        self.jsc = []
        self.vmax = []
        self.isc = []
        self.imax = []
        self.pmax = []
        self.fill_factor = []
        self.efficiency = []
        self.remainingFactorsColumnIndices = []

        for i, column in enumerate(self.qual_data.columns):
            if ('Voc' in column):
                if self.data[0,i] and ~np.isnan(self.data[0,i]):
                    self.voc = self.data[:, [0, 1, i]]
            elif ('Jsc' in column):
                if self.data[0,i] and ~np.isnan(self.data[0,i]):
                    self.jsc = self.data[:, [0, 1, i]]
            elif ('Vmax' in column):
                if self.data[0,i] and ~np.isnan(self.data[0,i]):
                    self.vmax = self.data[:, [0, 1, i]]
            elif ('Isc' in column):
                if self.data[0,i] and ~np.isnan(self.data[0,i]):
                    self.isc = self.data[:, [0, 1, i]]
            elif ('Imax' in column):
                if self.data[0,i] and ~np.isnan(self.data[0,i]):
                    self.imax = self.data[:, [0, 1, i]]
            elif ('Pmax' in column):
                if self.data[0,i] and ~np.isnan(self.data[0,i]):
                    self.pmax = self.data[:, [0, 1, i]]
                    # self.Efficiency = self.Pmax
            elif (column in ['FillFactor', 'FF', 'Fill Factor']):
                if self.data[0,i] and ~np.isnan(self.data[0,i]):
                    self.fill_factor = self.data[:, [0, 1, i]]
                    self.remaining_factors_available[i - 2] = u'FillFactor'  # minus 2 is because remaining Factor titles start after particle energy and fluence columns
            elif (column in ['Efficiency', 'eff', 'Eff', 'EFF']):
                if self.data[0,i] and ~np.isnan(self.data[0,i]):
                    self.efficiency = self.data[:, [0, 1, i]]
                    # self.Pmax = self.Efficiency
                    self.remaining_factors_available[i - 2] = u'Efficiency'

    def group_by_energy(self, param):
        data = getattr(self, param)
        indicesOfDuplicates = []
        energies = np.unique(data[:,0])
        for energy in energies:
            indicesOfDuplicates.append(list(np.where(data[:,0] == energy)[0]))
        grouped_by_energy = []
        for index in indicesOfDuplicates:
            grouped_by_energy.append(data[min(index):max(index) + 1])
        return grouped_by_energy
