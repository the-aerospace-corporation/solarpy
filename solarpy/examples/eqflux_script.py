import solarpy as sy
import pandas as pd

# Example Ae8/Ap8 data file  20,200km at 55 deg inclination
gps_10year_filename = r'C:\Users\dw28596\Documents\Data\Radiation Degradation Modeling\Radiation Spectra\GPS\GPS 10year orbit.xlsx'
gps_10year_data = pd.ExcelFile(gps_10year_filename)

differential_proton_spectrum = gps_10year_data.parse('Differential Protons').values[:, [0, 2]]
differential_electron_spectrum = gps_10year_data.parse('Differential Electrons').values[:, [0, 2]]

integral_proton_spectrum = gps_10year_data.parse('Integral Protons').values[:, [0,2]]
integral_electron_spectrum = gps_10year_data.parse('Integral Electrons').values[:, [0,2]]



# Interpolating differential particle spectrum according to EQFLUX method
interpolated_differential_proton_environment_spectrum = sy.eqflux.interpolate_proton_spectrum(differential_proton_spectrum)
interpolated_differential_electron_environment_spectrum = sy.eqflux.interpolate_electron_spectrum(differential_electron_spectrum)

# Don will provide appropriate RDC curves depending on cell technology as well as the proton to electron conversion factor

electron_rdc_curve_parameter = sy.GaAs_Electron_RDC[:,[0,5]] # This is just the GaAs RDC for fun
proton_rdc_curve_parameter = sy.GaAs_Proton_RDC_Pmax_Voc[:,[0,5]] # This is just the GaAs RDC for fun
proton_to_electron_conversion_factor = 1000 # Conversion factor for Pmax and Voc of GaAs RDC

# Below I put two methods to calculate the 1 MeV equivalent fluence.  The first using the differential particle spectra.
# I believe this method is more accurate than the original JPL EQFLUX method. It has to do with how the interpolate the
# RDC and use the integral spectrum.  In a lot of ways the original EQFLUX method seems overly complicated and
# unnecessary.  But I included it because it is the standard.  The first method tha uses the differential particle
# spectra and trapezoidal integration is cleaner and makes more mathmetical sense to me.  Your choice...also long as you
#  chose the same method.  Also for the differential method you have to interpolate the particle spectra before hand.
# The intperlation functions I use are in the library and basically interpolate how it is described in the JPL Radiation
#  handbooks, which are linear interpolations of the loglog proton spectra and linearlog of the electron spectra.

# Calculates 1 MeV Electron Fluence using the differential particle spectrum.
ten_MeV_protons_differential = sy.eqflux.get_relative_MeV_fluence_trapezoidal_integration(proton_rdc_curve_parameter, interpolated_differential_proton_environment_spectrum)
one_MeV_electrons_differential = sy.eqflux.get_relative_MeV_fluence_trapezoidal_integration(electron_rdc_curve_parameter, interpolated_differential_electron_environment_spectrum)
total_one_MeV_electrons_differential = one_MeV_electrons_differential + (ten_MeV_protons_differential * proton_to_electron_conversion_factor)

# Calculates the 1 MeV Electron Fluence using the exact method as in EQFLUX
ten_MeV_protons_integral = sy.eqflux.compute_equivalent_fluence_for_protons(proton_rdc_curve_parameter, integral_proton_spectrum)
one_MeV_electrons_integral = sy.eqflux.compute_equivalent_fluence_for_electrons(electron_rdc_curve_parameter, integral_electron_spectrum)
total_one_MeV_electrons_integral = one_MeV_electrons_integral + (ten_MeV_protons_integral * proton_to_electron_conversion_factor)

print('Total 1MeV electron fluence using differential particle spectra')
print(total_one_MeV_electrons_differential)
print()
print('Total 1Mev electron fluence using integral particle spectra')
print(total_one_MeV_electrons_integral)

