import os
import pandas as pd

class ImportEnvironmentSpectra(object):
    def __init__(self, environment_spectra_xlsx):
        """
        testing
        Args:
            environment_spectra_xlsx:
        """
        self.filename = os.path.basename(environment_spectra_xlsx)
        self.qualDataFrame = pd.ExcelFile(environment_spectra_xlsx)
        self.particle_environment_type = self.qualDataFrame.sheet_names
        self.electrons_integral = []
        self.protons_integral = []
        self.electrons_differential = []
        self.protons_differential = []
        self.solar_protons_integral = []
        self.solar_protons_differential = []


        for env in self.qualDataFrame.sheet_names:
            if env in ['Electrons Integral']:
                self.electrons_integral = self.qualDataFrame.parse(env).values.astype('float')
            elif env in ['Protons Integral']:
                self.protons_integral = self.qualDataFrame.parse(env).values.astype('float')
            elif env in ['Electrons Differential']:
                self.electrons_differential = self.qualDataFrame.parse(env).values.astype('float')
            elif env in ['Protons Differential']:
                self.protons_differential = self.qualDataFrame.parse(env).values.astype('float')
            elif env in ['Solar Protons Integral']:
                self.solar_protons_integral = self.qualDataFrame.parse(env).values.astype('float')
            elif env in ['Solar Protons Differential']:
                self.solar_protons_differential = self.qualDataFrame.parse(env).values.astype('float')
