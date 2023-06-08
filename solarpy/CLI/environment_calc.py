import argparse
import pandas as pd
import json
import os
import numpy as np
from solarpy import data
import solarpy.ddd.ddd_degradation_functions as ddd_eq
import solarpy.degradation_equations as eq
from solarpy import file_imports as sy_imports
import solarpy.eqflux.relativeMeV_fluence as rel_MeV

params = ['pmax', 'voc', 'vmax', 'isc', 'imax', 'fillfactor', 'efficiency']

def main():
    parser = argparse.ArgumentParser(prog='Environment Calculator')
    parser.add_argument('-p', '--proton_env', help='Proton environment spectra', metavar='')
    parser.add_argument('-e', '--electron_env', help='Electron environment spectra', metavar='')
    parser.add_argument('-et', '--environment_type', help='Whether environment type is integral or differential (differential prefered)', metavar='int or diff', default='diff')
    parser.add_argument('-param', help='Electrical Parameter', choices=params)
    parser.add_argument('-prdc', '--proton_rdc', help='Proton RDC curve to be used for EQFLUX calculations', metavar='')
    parser.add_argument('-erdc', '--electron_rdc', help='Electron RDC curve to be used for EQFLUX calculations', metavar='')
    parser.add_argument('-pe', '--p_to_e_conversion', help='Conversion factor to convert protons to electrons', metavar='')
    parser.add_argument('-ep', '--e_to_p_conversion', help='Conversion factor to convert electrons to protons', metavar='')
    parser.add_argument('-pd', '--proton_ddd', help='Proton DDD fit parameters', default=None, nargs=2, type=float, metavar='C Phi')
    parser.add_argument('-ed', '--electron_ddd', help='Proton DDD fit parameters', default=None, nargs=2, type=float, metavar='C Phi')

    parser.add_argument('-pd', '--proton_ddd', help='Starting guess to fit C and phi or to calculate remaining factor', default=None, nargs='*', type=float, metavar='C Phi A')
    parser.add_argument('-dfp', '--double_fit_param', help='Starting guess to fit C1, phi1, C2, phi2 or to calculate remaining factor', default=None, nargs='*', type=float, metavar='C1 phi1 C2 phi2 A')
    args = parser.parse_args()

    rel_MeV_total = []

    if args.proton_env and args.proton_rdc:
        env = sy_imports.get_dataframe_from_input_file(args.proton_env)
        if args.environment_type is 'diff' or (env['type'][0] == 'differential'):
            proton_env = env[['energy_mev', 'fluence']].values
            intp_diff_spectrum = rel_MeV.interpolate_proton_spectrum(proton_env)
            rdc = sy_imports.get_dataframe_from_input_file(args.proton_rdc).values
            rel_MeV_total = rel_MeV.get_relative_MeV_fluence_trapezoidal_integration(rdc, intp_diff_spectrum)

    if args.electron_env and args.electron_rdc:
        env = sy_imports.get_dataframe_from_input_file(args.electron_env)
        if args.environment_type is 'diff' or (env['type'][0] == 'differential'):
            electron_env = env[['energy_mev', 'fluence']].values
            intp_diff_spectrum = rel_MeV.interpolate_electron_spectrum(electron_env)
            rdc = sy_imports.get_dataframe_from_input_file(args.electron_rdc).values
            rel_MeV_total = rel_MeV.get_relative_MeV_fluence_trapezoidal_integration(rdc, intp_diff_spectrum)

    if args.proton_env and args.electron_rdc and args.proton_rdc and args.electron_rdc and (args.p_to_e_conversion or args.e_to_p_conversion):
        env_p = sy_imports.get_dataframe_from_input_file(args.proton_env)
        env_e = sy_imports.get_dataframe_from_input_file(args.electron_env)

        if args.environment_type is 'diff' or (env_p['type'][0] == 'differential'):
            electron_env = env_e[['energy_mev', 'fluence']].values
            proton_env = env_p[['energy_mev', 'fluence']].values

            intp_diff_spectrum_e = rel_MeV.interpolate_electron_spectrum(electron_env)
            intp_diff_spectrum_p = rel_MeV.interpolate_proton_spectrum(proton_env)

            rdc_e = sy_imports.get_dataframe_from_input_file(args.electron_rdc).values
            rdc_p = sy_imports.get_dataframe_from_input_file(args.proton_rdc).values

            electron_rel_MeV_total = rel_MeV.get_relative_MeV_fluence_trapezoidal_integration(rdc_e, intp_diff_spectrum_e)
            proton_rel_MeV_total = rel_MeV.get_relative_MeV_fluence_trapezoidal_integration(rdc_p, intp_diff_spectrum_p)

            if args.p_to_e_conversion:
                rel_MeV_total = electron_rel_MeV_total + (proton_rel_MeV_total*args.p_to_e_conversion)

            elif args.e_to_p_conversion:
                rel_MeV_total = proton_rel_MeV_total + (electron_rel_MeV_total*args.e_to_p_conversion)

    # if args.proton_env and args.proton_ddd:

                