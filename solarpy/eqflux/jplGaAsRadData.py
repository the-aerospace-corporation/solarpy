import pandas as pd
import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import solarpy as sy

class jplGaAsRadData:
    def __init__(self, JPL_GaAs_Radiation_Handbook_Data_xlsx):
        self.file = JPL_GaAs_Radiation_Handbook_Data_xlsx
        self.allData = pd.ExcelFile(self.file, engine='openpyxl')
        self.rawProtons = self.allData.parse('Protons (extracted)').values[1:].astype('float')
        self.rawElectrons = self.allData.parse('Electrons (extracted)').values[1:].astype('float')
        self.protonEnergies = np.array([0.05, 0.1, 0.2, 0.3, 0.5, 1, 3, 9.5])
        self.electronEnergies = np.array([0.6, 1, 2.4, 12])
        self.rawElectronByEnergy = self._groupedByEnergy_Electrons()
        self.rawProtonByEnergy = self._groupedByEnergy_Protons()
        self.protonQualData = np.vstack(self.getProtonQualDataGroupedByenergy())
        self.electronQualData = np.vstack(self.getElectronQualDataGroupedByEnergy())

    def getProtonQualDataGroupedByenergy(self):

        proton50keV = np.array([3e9, 5e9, 1e10, 3e10, 1e11, 1.8e11, 2.8e11]) #np.array([1e10, 2e10, 4e10, 6e10, 8e10, 1e11])
        proton100keV = np.array([5e9, 1e10, 5e10, 1e11, 1.4e11, 2.7e11])#np.array([1e10, 2e10, 4e10, 6e10, 8e10, 1e11])
        proton200keV = np.array([3e9, 5e9, 1e10, 3e10, 5e10, 1e11, 2e11, 2.7e11])#np.array([1e10, 2e10, 4e10, 6e10, 8e10, 1e11, 1.5e11, 2e11])
        proton300keV =  np.array([3e9, 5e9, 1e10, 3e10, 5e10, 1e11, 2.2e11, 2.7e11])#np.array([1e10, 2e10, 4e10, 6e10, 8e10, 1e11, 1.5e11, 2e11])
        proton500keV =  np.array([3e9, 5e9, 1e10, 3e10, 5e10, 1e11,9.4e11])#np.array([1e10, 2e10, 4e10, 6e10, 8e10, 1e11, 1.5e11, 2e11, 4e11, 8e11])
        proton1MeV =  np.array([1e10, 2e10, 5e10, 1e11, 3e11, 9.3e11, 9.5e11])#np.array([1e10, 2e10, 4e10, 6e10, 8e10, 1e11, 1.5e11, 2e11, 4e11, 8e11])
        proton3MeV = np.array([5e10, 1e11, 5e11, 1e12, 2e12, 2.8e12, 2.9e12])#np.array([1e10, 2e10, 4e10, 6e10, 8e10, 1e11, 1.5e11, 2e11, 4e11, 8e11, 1e12, 2e12])
        proton9_5MeV = np.array([5e10, 1e11, 5e11, 1e12, 2e12, 2.6e12, 9.7e12])#np.array([1e10, 2e10, 4e10, 6e10, 8e10, 1e11, 1.5e11, 2e11, 4e11, 8e11, 1e12, 2e12])

        dataTable = []
        for i, energy in enumerate(self.protonEnergies):
            if  energy == 0.05:
                fluence = proton50keV
            elif energy == 0.1:
                fluence = proton100keV
            elif energy == 0.2:
                fluence = proton200keV
            elif energy == 0.3:
                fluence = proton300keV
            elif energy == 0.5:
                fluence = proton500keV
            elif energy == 1:
                fluence = proton1MeV
            elif energy == 3:
                fluence = proton3MeV
            elif energy == 9.5:
                fluence = proton9_5MeV

            # print energy
            # print fluence
            # print self.rawProtonByEnergy[i]
            fluence = fluence[fluence>=np.min(self.rawProtonByEnergy[i][:,0])]
            fluence = fluence[fluence<=np.max(self.rawProtonByEnergy[i][:,0])]
            particleEnergyArray = np.array([energy] * (len(fluence) ))
            data = self.intpByFluence(fluence,self.rawProtonByEnergy[i])
            data = np.column_stack((particleEnergyArray.T,data[:,0], data[:,1]))
            dataTable.append(data)
        return dataTable

    def getElectronQualDataGroupedByEnergy(self):

        electron600keV = np.array([2e13, 2e14, 1e15, 3e15, 1e16])
        electron1MeV = np.array([5e13, 1e14, 2e14, 5e14, 1e15, 3e15, 1e16])
        electron2_4MeV = np.array([5e12, 5e13, 2e14, 5e14, 1e15, 3.5e15])
        electron12MeV = np.array([2e12, 1e13, 5e13, 1e14, 3e14, 5e14, 1.2e15])

        # electron600keV = np.array([2e14, 1e15, 3e15, 1e16])
        # electron1MeV = np.array([2e14, 5e14, 1e15, 3e15, 1e16])
        # electron2_4MeV = np.array([2e14, 5e14, 1e15, 3.5e15])
        # electron12MeV = np.array([2e14, 4e14, 8e14, 2e15])


        dataTable = []
        for i, energy in enumerate(self.electronEnergies):
            if  energy == 0.6:
                fluence = electron600keV
            elif energy == 1:
                fluence = electron1MeV
            elif energy == 2.4:
                fluence = electron2_4MeV
            elif energy == 12:
                fluence = electron12MeV
            fluence = fluence[fluence >= np.min(self.rawElectronByEnergy[i][:, 0])]
            fluence = fluence[fluence <= np.max(self.rawElectronByEnergy[i][:, 0])]
            particleEnergyArray = np.array([energy] * (len(fluence) ))
            data = self.intpByFluence(fluence,self.rawElectronByEnergy[i])
            data = np.column_stack((particleEnergyArray.T,data[:,0], data[:,1]))
            dataTable.append(data)
        return dataTable

    def intpByFluence(self, fluence, fluenceVSremainingFactor):
        newFluenceVSrf = np.zeros((len(fluence),2))
        newFluenceVSrf[:,0] = fluence
        intpRF = sp.interpolate.interp1d(fluenceVSremainingFactor[:,0], fluenceVSremainingFactor[:,1])
        newFluenceVSrf[:,1] = intpRF(fluence)
        # newFluenceVSrf = np.vstack((np.array([0,1]),newFluenceVSrf))
        return newFluenceVSrf

    def _groupedByEnergy_Protons(self):
        groupedByEnergy = []
        j=0
        for i in range(len(self.protonEnergies)):
            data = self.rawProtons[:,[j,j+1]]
            data = data[~np.isnan(data[:,0])]
            # data = pd.to_numeric(data)
            # print data
            groupedByEnergy.append(data)
            j = j+2
        groupedByEnergy = np.array(groupedByEnergy, dtype=object)
        return groupedByEnergy

    def _groupedByEnergy_Electrons(self):
        groupedByEnergy = []
        j=0
        for i in range(len(self.electronEnergies)):
            data = self.rawElectrons[:,[j,j+1]]
            data = data[~np.isnan(data[:,0])]
            # data = pd.to_numeric(data)
            # print data
            groupedByEnergy.append(data)
            j = j+2
        groupedByEnergy = np.array(groupedByEnergy, dtype=object)
        return groupedByEnergy

    def plotJPLGaAsData(self):
        electronPlot = []
        electronCubicSpline = self.rawElectronByEnergy
        protonCubicSpline = self.rawProtonByEnergy
        colors = ['blue', 'g', 'r', 'c', 'purple', 'yellow', 'black', 'orange']
        for i, electrons in enumerate(self.getElectronQualDataGroupedByEnergy()):
            l = str(electrons[0,0])+r' MeV'
            electronPlot.append(plt.semilogx(electrons[:,1], electrons[:,2], 'o', mec=None, color=colors[i],label= l))
            plt.semilogx(electronCubicSpline[i][:,0], electronCubicSpline[i][:,1], lw = 2, color=colors[i])
        legend1 = plt.legend(handles=[label[0] for label in electronPlot], loc=4, title="Electrons", frameon=False, numpoints=1, prop={'size':12, 'family':'Arial'})
        protonPlot = []
        for i, protons in enumerate(self.getProtonQualDataGroupedByenergy()):
            l = str(protons[0,0])+' MeV'
            protonPlot.append(plt.semilogx(protons[:,1], protons[:,2], 's', mec=None, color=colors[i], label= l))
            plt.semilogx(protonCubicSpline[i][:,0], protonCubicSpline[i][:,1], lw = 2, color=colors[i])
        legend2 = plt.legend(handles=[label[0] for label in protonPlot], loc=3, title='Protons',frameon=False, numpoints=1, prop={'size':12, 'family':'Arial'})
        plt.gca().add_artist(legend1)
        # sy.plotSetup(xlabel=r'Fluence (cm$^2$)', ylabel=r'P$_{max}$/P$_0$',y_limits=[0.1, 1.05])
        plt.xscale('symlog', linthreshx = 1.1e9)
        plt.xlim([0,1.2e16])