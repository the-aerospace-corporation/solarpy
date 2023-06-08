import solarpy.ddd.ddd_degradation_functions as ddd_eq
import solarpy.degradation_equations as d_eq
from solarpy import data
import solarpy.eqflux.relativeMeV_fluence as oneMeV

class ddd_aero(object):
    def __init__(self, electron_qual_data, proton_qual_data, proton_slowed_down_differential_spectrum, electron_slowed_down_differential_spectrum):
        self.electron_qual_data = electron_qual_data
        self.proton_qual_data = proton_qual_data
        self.fully_penetrating_qual_data = self.proton_qual_data[self.proton_qual_data[:,0] >= 1]
        self.proton_energy_vs_ddd = ddd_eq.get_ddd_vs_remaining_factor(self.fully_penetrating_qual_data[:,0],self.fully_penetrating_qual_data[:,1],self.fully_penetrating_qual_data[:,2],data.proton_NIEL_SCREAM)
        self.n = ddd_eq.fit_nValue_effectiveDDD(self.electron_qual_data[:,0], self.electron_qual_data[:,1], self.electron_qual_data[:,2], data.electron_NIEL_SCREAM, energy_to_normalize=1)
        self.electron_energy_vs_ddd = ddd_eq.get_ddd_vs_remaining_factor(self.electron_qual_data[:,0], self.electron_qual_data[:,1], self.electron_qual_data[:,2], data.electron_NIEL_SCREAM, energy_to_normalize=1, n=self.n)
        self.electron_reference_object = d_eq.singleDegradationEquation(self.electron_energy_vs_ddd)
        self.proton_reference_object = d_eq.singleDegradationEquation(self.proton_energy_vs_ddd)
        self.proton_ddd_spectrum = ddd_eq.get_energy_vs_ddd(oneMeV.interpolate_proton_spectrum(proton_slowed_down_differential_spectrum), data.proton_NIEL_SCREAM)
        self.electron_ddd_spectrum = ddd_eq.get_energy_vs_ddd(oneMeV.interpolate_electron_spectrum(electron_slowed_down_differential_spectrum), data.electron_NIEL_SCREAM)
        self.proton_total_ddd = ddd_eq.getTotalDDD(self.proton_ddd_spectrum)
        self.electron_total_ddd = ddd_eq.getTotalDDD(self.electron_ddd_spectrum)
        self.totalDDD_from_spectra = self.proton_total_ddd + self.convertDDDelectronsToDDDprotons(self.electron_total_ddd)
        self.remaining_factor = self.proton_reference_object.getRemainingFactor(self.totalDDD_from_spectra)

    def get_reference_ddd_object(self, energy_vs_ddd):
        return d_eq.singleDegradationEquation(energy_vs_ddd)

    def convertDDDprotonsToDDDelectrons(self, D_p):
        self.electron_reference_object.fit()
        self.proton_reference_object.fit()

        electronFit = self.electron_reference_object.coefficients
        protonFit = self.proton_reference_object.coefficients
        C_e = electronFit[0]
        C_p = protonFit[0]
        D_ex = electronFit[1]
        D_px = protonFit[1]
        D_pTOe = ddd_eq.convertDDDprotonsToDDDelectrons(D_px,D_ex, C_e, C_p, D_p)
        return D_pTOe

    def convertDDDelectronsToDDDprotons(self, D_e):
        self.electron_reference_object.fit()
        self.proton_reference_object.fit()

        electronFit = self.electron_reference_object.coefficients
        protonFit = self.proton_reference_object.coefficients
        C_e = electronFit[0]
        C_p = protonFit[0]
        D_ex = electronFit[1]
        D_px = protonFit[1]
        D_eTOp =  ddd_eq.convertDDDelectronsToDDDprotons(D_px, D_ex, C_e, C_p, D_e)
        return D_eTOp



