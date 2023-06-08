import os
import numpy as np


class ImportTrappedRadiationSpenvis(object):
    def __init__(self, spenvis_tri_txt):
        self.file_name = os.path.basename(spenvis_tri_txt)
        # self.proton_spectra = []
        # self.electron_spectra = []
        self.eof = 1

        with open(spenvis_tri_txt, 'rU') as f:
            while self.eof:
                self.get_data(f)


    def get_data(self, file_object):
        trapped_proton_model = False
        trapped_electron_model = False
        mission_duration_days = []
        data = []
        for j in file_object:
            m = [a.strip().strip('\'') for a in j.rstrip().split(':')]
            n = [a.strip().strip('\'') for a in j.rstrip().split(',')]
            o = m + n
            if 'Trapped proton model' in o:
                trapped_proton_model = True
            elif 'Trapped electron model' in o:
                trapped_electron_model = True
            elif 'MIS_DUR' in o:
                mission_duration_days = float(n[2])
            elif 'Differential Flux' in o:
                break
            elif 'End of File' in o:
                self.eof = 0
                break

        if self.eof:
            for j in file_object:
                m = [k.strip().strip('\'') for k in j.rstrip().split(',')]
                if 'End of Block' in m:
                    break
                elif 'End of File' in m:
                    self.eof = 0
                    break
                else:
                    data.append([float(n) for n in m])
            if trapped_proton_model:
                self.proton_spectra = TrappedRadiationSpectra(np.asarray(data))
                self.proton_spectra.get_fluence(mission_duration_days)
                self.proton_spectra.particle = 'p'

            elif trapped_electron_model:
                self.electron_spectra = TrappedRadiationSpectra(np.asarray(data))
                self.electron_spectra.get_fluence(mission_duration_days)
                self.electron_spectra.particle = 'e'


class TrappedRadiationSpectra(object):
    def __init__(self, rad_data):
        self.particle = []
        self.mission_duration_days = []
        self.energy = rad_data[0] #MeV
        self.integral_flux = rad_data[:,[0,1]] # #cm^2/s
        self.differential_flux = rad_data[:,[0,2]] # #cm^2/s/MeV
        self.particle_spectra = rad_data
        self.integral_spectra = [] # Energy (MeV), Flux, Fluence (mission)
        self.differential_spectra = [] # Energy (MeV), Flux, Fluence (mission)

    def get_fluence(self, mission_duration_days=None):
        if mission_duration_days is not None:
            self.mission_duration_days = mission_duration_days

        if self.mission_duration_days:
            self.integral_spectra = np.vstack((self.integral_flux.T, self.integral_flux[:,1]*self.mission_duration_days * 24 * 60 * 60)).T
            self.differential_spectra = np.vstack((self.differential_flux.T, self.differential_flux[:,1]*self.mission_duration_days * 24 * 60 * 60)).T

