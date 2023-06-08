import os
import types

import numpy as np


class ImportMulassisSpenvisOutputFile(object): #TODO: make sure you make a bucket for electron and proton
    def __init__(self, spenvis_mlo_txt):
        self.file_name = os.path.basename(spenvis_mlo_txt)

        # self.electron_omnidirectional_fluence = []
        # self.electron_omnidirectional_fluence_0_90 = []
        # self.electron_omnidirectional_fluence_90_180 = []
        #
        # self.proton_omnidirectional_fluence = []
        # self.proton_omnidirectional_fluence_0_90 = []
        # self.proton_omnidirectional_fluence_90_180 = []
        self.eof = 1
        with open(spenvis_mlo_txt, 'rU') as f:
            while self.eof:
                data = self.get_data(f)


    def get_data(self, file_object):
        omni_flag = False
        omni_0_90_flag = False
        omni_90_180_flag = False
        proton_flag = False
        electron_flag = False
        flu_amn = False
        flu_amx = False

        for j in file_object:
            k = j.rstrip().split(',')
            m = [a.strip().strip('\'') for a in k]
            if any(terms in m for terms in ['OMNIDIRECTIONAL FLUENCE ANALYSIS', 'OMNIDIRECTIONAL FLUENCE ANALYSIS AS A FUNCTION OF ANGLE']):
                omni_flag = True
            elif 'proton' in m or 'p+' in m:
                proton_flag = True
            elif 'electron' in m or 'e-' in m:
                electron_flag = True
            elif 'FLU_AMN' in m:
                flu_amn = float(m[2])
            elif 'FLU_AMX' in m:
                flu_amx = float(m[2])
            elif 'Error' in m:
                break
            elif 'End of File' in m:
                self.eof = 0
                break

        if omni_flag:
            if flu_amn == 0 and flu_amx == 90:
                omni_flag = False
                omni_0_90_flag = True
            elif flu_amn == 90 and flu_amx == 180:
                omni_flag = False
                omni_90_180_flag = True
        data = []

        if self.eof:
            for j in file_object:
                k = j.rstrip().split(',')
                m = [a.strip().strip('\'') for a in k]

                if 'End of Block' in m:
                    break
                elif 'End of File' in m:
                    self.eof = 0
                    break
                else:
                    data.append([float(n) for n in m])

            spectra = mulassis_spectra(np.asarray(data))
            spectra.flu_amn = flu_amn
            spectra.flu_amx = flu_amx

            if proton_flag:
                spectra.particle = 'p'
                if omni_flag:
                    self.proton_omnidirectional_fluence = spectra
                elif omni_0_90_flag:
                    self.proton_omnidirectional_fluence_0_90 = spectra
                elif omni_90_180_flag:
                    self.proton_omnidirectional_fluence_90_180 = spectra
            elif electron_flag:
                spectra.particle = 'e'
                if omni_flag:
                    self.electron_omnidirectional_fluence = spectra
                elif omni_0_90_flag:
                    self.electron_omnidirectional_fluence_0_90 = spectra
                elif omni_90_180_flag:
                    self.electron_omnidirectional_fluence_90_180 = spectra



class get_flags():
    def __init__(self, file_object):

        self.omni_flag = False
        self.omni_0_90_flag = False
        self.omni_90_180_flag = False
        self.proton_flag = False
        self.electron_flag = False
        self.flu_amn = False
        self.flu_amx = False
        self.eof = 1

        for j in file_object:
            k = j.rstrip().split(',')
            m = [a.strip().strip('\'') for a in k]
            if any(terms in m for terms in ['OMNIDIRECTIONAL FLUENCE ANALYSIS', 'OMNIDIRECTIONAL FLUENCE ANALYSIS AS A FUNCTION OF ANGLE']):
                self.omni_flag = True
            elif 'proton' in m or 'p+' in m:
                self.proton_flag = True
            elif 'electron' in m or 'e-' in m:
                self.electron_flag = True
            elif 'FLU_AMN' in m:
                self.flu_amn = float(m[2])
            elif 'FLU_AMX' in m:
                self.flu_amx = float(m[2])
            elif 'Error' in m:
                break
            elif 'End of File' in m:
                self.eof = 0
                break

        if self.omni_flag:
            if self.flu_amn == 0 and self.flu_amx == 90:
                self.omni_flag = False
                self.omni_0_90_flag = True
            elif self.flu_amn == 90 and self.flu_amx == 180:
                self.omni_flag = False
                self.omni_90_180_flag = True

class mulassis_spectra:
    def __init__(self, mulassis_data):
        mulassis_data = mulassis_data[mulassis_data[:,2]>0]
        self.particle = []
        self.flu_amn = []
        self.flu_amx = []
        self.elo = mulassis_data[:,0] # keV
        self.eup = mulassis_data[:,1] # keV
        self.emean = mulassis_data[:,2] # keV
        self.value = mulassis_data[:,3]
        self.error = mulassis_data[:,4]
        self.slowed_down_spectra = mulassis_data[:,[2,3]]
        self.data = mulassis_data



