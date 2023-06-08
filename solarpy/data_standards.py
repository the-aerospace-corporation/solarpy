import pandas as pd
import numpy as np
from solarpy.eqflux.jplGaAsRadData import jplGaAsRadData
import pkg_resources

# GaAs RDC from The GaAs Radiation Handbook
# column 0 is particle energy, column 1 is normal incident RDC, column 2 is omnidirectional RDC for no shielding and columns 3-9 are RDCS for CMG shieldings of thickness 1 mil, 3 mil, 6 mll, 12 mil, 20 mil, 30 mil, and 60 mil

GaAs_Electron_RDC = pd.ExcelFile(pkg_resources.resource_filename('solarpy', 'data/GaAs Damage Coefficiencts (JPL GaAs Radiation Handbook).xlsx'), engine='openpyxl').parse('Electron Damage Coeffecients').values
GaAs_Proton_RDC_Pmax_Voc = pd.ExcelFile(pkg_resources.resource_filename('solarpy', 'data/GaAs Damage Coefficiencts (JPL GaAs Radiation Handbook).xlsx'), engine='openpyxl').parse('Proton Pmax and Voc Damage Coef').values
GaAs_Proton_RDC_Isc = pd.ExcelFile(pkg_resources.resource_filename('solarpy', 'data/GaAs Damage Coefficiencts (JPL GaAs Radiation Handbook).xlsx'), engine='openpyxl').parse('Proton Isc Damage Coeffecients').values

# Silicon RDC from The Solar Cell Radiation Handbook
# column 0 is particle energy, column 1 is omnidirectional RDC for no shielding and columns 2-8 are RDCS for CMG shieldings of thickness 1 mil, 3 mil, 6 mll, 12 mil, 20 mil, 30 mil, and 60 mil

Si_Electron_RDC = pd.ExcelFile(pkg_resources.resource_filename('solarpy', 'data/Si Damage Coefficiencts (JPL GaAs Radiation Handbook).xlsx'), engine='openpyxl').parse('Electron Damage Coeffecients').values
Si_Proton_RDC_Pmax_Voc = pd.ExcelFile(pkg_resources.resource_filename('solarpy', 'data/Si Damage Coefficiencts (JPL GaAs Radiation Handbook).xlsx'), engine='openpyxl').parse('Proton Pmax and Voc Damage Coef').values
Si_Proton_RDC_Isc = pd.ExcelFile(pkg_resources.resource_filename('solarpy', 'data/Si Damage Coefficiencts (JPL GaAs Radiation Handbook).xlsx'), engine='openpyxl').parse('Proton Isc Damage Coeffecients').values

#GaAs RDC from The GaAs Radiation Handbook calculated using Spenvis RDC calculator
# column 0 is particle energy, column 1 is normal incident RDC, column 2 is omnidirectional RDC for no shielding and columns 3-9 are RDCS for CMG shieldings of thickness 1 mil, 3 mil, 6 mll, 12 mil, 20 mil, 30 mil, and 60 mil

GaAs_Electron_RDC_Spenvis = pd.ExcelFile(pkg_resources.resource_filename('solarpy', 'data/GaAs Damage Coefficiencts (Spenvis).xlsx'), engine='openpyxl').parse('Electron Damage Coeffecients').values
GaAs_Proton_RDC_Pmax_Voc_Spenvis = pd.ExcelFile(pkg_resources.resource_filename('solarpy', 'data/GaAs Damage Coefficiencts (Spenvis).xlsx'), engine='openpyxl').parse('Proton Pmax Damage Coefficients').values
GaAs_Proton_RDC_Isc_Spenvis = pd.ExcelFile(pkg_resources.resource_filename('solarpy', 'data/GaAs Damage Coefficiencts (Spenvis).xlsx'), engine='openpyxl').parse('Proton Isc Damage Coefficients').values

#CMG Range Table
CMG_proton_ranges = pd.ExcelFile(pkg_resources.resource_filename('solarpy', 'data/JPL CMG Range Tables.xlsx'), engine='openpyxl').parse('Protons').values
CMG_electron_ranges = pd.ExcelFile(pkg_resources.resource_filename('solarpy', 'data/JPL CMG Range Tables.xlsx'), engine='openpyxl').parse('Electrons').values
FS_proton_ranges = pd.ExcelFile(pkg_resources.resource_filename('solarpy', 'data/JPL CMG Range Tables.xlsx'), engine='openpyxl').parse('Protons').values
FS_electron_ranges = pd.ExcelFile(pkg_resources.resource_filename('solarpy', 'data/JPL CMG Range Tables.xlsx'), engine='openpyxl').parse('Electrons').values
CMG_density = 2.554
FS_density = 2.2

#GaAs Range Table and Displacements
GaAs_proton_ranges = pd.ExcelFile(pkg_resources.resource_filename('solarpy', 'data/protons in GaAs Range Table.xlsx'), engine='openpyxl').parse('Proton Range Table GaAs').values
GaAs_total_displacements = pd.ExcelFile(pkg_resources.resource_filename('solarpy', 'data/protons in GaAs Range Table.xlsx'), engine='openpyxl').parse('TRIM Displacements').values

#GaAs Qual Data as extracted from xxxx reference
gaas_qual_data = jplGaAsRadData(pkg_resources.resource_filename('solarpy', 'data/JPL GaAs Radiation Handbook Power_High_Res.xlsx'))
gaas_qual_data_power = jplGaAsRadData(pkg_resources.resource_filename('solarpy', 'data/JPL GaAs Radiation Handbook Power_High_Res.xlsx'))
gaas_qual_data_isc = jplGaAsRadData(pkg_resources.resource_filename('solarpy', 'data/JPL GaAs Radiation Handbook Isc_High_Res.xlsx'))
gaas_qual_data_voc = jplGaAsRadData(pkg_resources.resource_filename('solarpy', 'data/JPL GaAs Radiation Handbook Voc_High_Res.xlsx'))

#NIEL curves for GaAs from SCREAM
proton_NIEL_SCREAM = np.loadtxt(pkg_resources.resource_filename('solarpy', 'data/GaAs_Scream_Protons.txt'),skiprows=1)
electron_NIEL_SCREAM = np.loadtxt(pkg_resources.resource_filename('solarpy', 'data/GaAs_Scream_Electrons.txt'),skiprows=1)
