import numpy as np
import scipy as sp
import solarpy.eqflux.relativeMeV_fluence as rel_MeV


class EQFLUX_electrons(object):
    def __init__(self, electron_energy_vs_rdc, electron_differential_spectrum):
        self.electron_energy_vs_rdc = electron_energy_vs_rdc #typically relative to 1 MeV
        self.differential_electron_spectrum = electron_differential_spectrum
        self.interpolated_electron_spectrum = rel_MeV.interpolate_electron_spectrum(electron_differential_spectrum)

    def rel_MeV_Fluence(self):
        oneMeVFluence = rel_MeV.get_relative_MeV_fluence_trapezoidal_integration(self.electron_energy_vs_rdc, self.interpolated_electron_spectrum)
        return oneMeVFluence

    def get_energy_vs_rel_MeV_fluence(self):
        energy_vs_rel_MeV_fluence = rel_MeV.get_energy_vs_relative_fluence(self.electron_energy_vs_rdc, self.interpolated_electron_spectrum)
        return energy_vs_rel_MeV_fluence

    def get_cumulative_rel_MeV_fluence(self):
        cumulative_rel_MeV_fluence = rel_MeV.get_cumulative_relative_fluence(self.electron_energy_vs_rdc, self.interpolated_electron_spectrum)
        return cumulative_rel_MeV_fluence
