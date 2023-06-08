import numpy as np
import scipy as sp
import solarpy.eqflux.relativeMeV_fluence as rel_MeV


class EQFLUX_protons(object):
    def __init__(self, proton_energy_vs_rdc, proton_differential_spectrum):
        self.proton_energy_vs_rdc = proton_energy_vs_rdc #typically relative to 10 MeV
        self.differential_proton_spectrum = proton_differential_spectrum
        self.interpolated_proton_spectrum = rel_MeV.interpolate_proton_spectrum(proton_differential_spectrum)

    def ten_MeV_Fluence(self):
        oneMeVFluence = rel_MeV.get_relative_MeV_fluence_trapezoidal_integration(self.proton_energy_vs_rdc, self.interpolated_proton_spectrum)
        return oneMeVFluence

    def get_energy_vs_10MeV_fluence(self):
        energy_vs_1MeV_fluence = rel_MeV.get_energy_vs_relative_fluence(self.proton_energy_vs_rdc, self.interpolated_proton_spectrum)
        return energy_vs_1MeV_fluence

    def get_cumulative_10MeV_fluence(self):
        cumulative_1MeV_fluence = rel_MeV.get_cumulative_relative_fluence(self.proton_energy_vs_rdc, self.interpolated_proton_spectrum)
        return cumulative_1MeV_fluence
