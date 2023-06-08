import argparse
import pandas as pd
import json
import os
from solarpy import data
import solarpy.ddd.ddd_degradation_functions as ddd_eq
import solarpy.degradation_equations as eq
from solarpy import file_imports as sy_imports

electron_rdc_vars = ['electron qual data', 'electrical parameter', 'critical remaining factor',
                     'energy to normalize rdc', 'shield thickness in cm', 'density of shield', 'shield range table']
proton_rdc_vars = ['proton qual data', 'electrical parameter', 'critical remaining factor', 'energy to normalize rdc',
                   'shield thickness in cm', 'solar cell thickness' 'density of shield', 'shield range table',
                   'solar cell range table', 'solar cell displacement table']
params = ['pmax', 'voc', 'vmax', 'isc', 'imax', 'fillfactor', 'efficiency']


def main():
    parser = argparse.ArgumentParser(prog='DDD Calculator')
    parser.add_argument('-p', '--particle',
                        help='Calculate RDC for particle type: protons (p), electrons (e), both (b)',
                        choices=['p', 'e', 'b'])
    parser.add_argument('-param', help='Electrical Parameter', choices=params)
    parser.add_argument('-f', '--fit_type',
                        help='Select fit type: Single Degradation Equation (s), Double Degradation Equation(d)',
                        choices=['single', 'double'], default='s')
    parser.add_argument('-sfp', '--single_fit_param', help='Starting guess to fit C and phi or to calculate remaining factor', default=None, nargs=2, type=float, metavar='C Phi A')
    parser.add_argument('-dfp', '--double_fit_param', help='Starting guess to fit C1, phi1, C2, phi2 or to calculate remaining factor', default=None, nargs=4, type=float, metavar='C1 phi1 C2 phi2 A')
    parser.add_argument('-r', '--rad', help='Radiation Data in either xlsx or txt', metavar='')
    parser.add_argument('-c', '--cut_off', help='Proton energy above which particles are fully penetrating', default=1, metavar='')
    parser.add_argument('-norm_e', help='Electron energy to normalize electron radiation data', type=float, default=1,
                        metavar='')
    parser.add_argument('-norm_p', help='Proton energy to normalize proton radiation data', type=float, default=10,
                        metavar='')
    parser.add_argument('-rf', '--remaining_factor', help='Input desired DDD to get remaining factor', default = None, metavar='', type=float)
    parser.add_argument('-st', help='Shield thickness in cm', type=float, default=0, metavar='')
    parser.add_argument('-sd', help='Shield density in g/cm2', type=float, default=2.554, metavar='')
    parser.add_argument('-srt', help='Shield range table xls', metavar='')
    parser.add_argument('-od', '--output_loc', help='Where to put output files', metavar='')
    parser.add_argument('-o', '--output_xlsx', help='Output Table. Give filename', metavar='')
    parser.add_argument('-oj', '--output_json', help='Output json. Give filename', metavar='')
    parser.add_argument('-oa', '--output_all', help='Output json. Give filename', metavar='')

    args = parser.parse_args()

    #initialize data lists
    ddd_vs_rf = []
    qual_data_frame = []
    coefficients = []
    fits = []
    rf=[]

    ### takes in input data as txt, json, xls and starts to calculate ddd vs rf and fit it to user specified single or double exponential fits
    if args.rad:
        qual_data_frame = sy_imports.get_dataframe_from_input_file(args.rad, type='qual')
        # print(qual_data_frame)
        if args.particle is 'p':
            qual_data_frame = qual_data_frame[qual_data_frame.energy_mev >= args.cut_off]
            ddd_vs_rf = ddd_eq.get_ddd_vs_remaining_factor(particle_energy=qual_data_frame['energy_mev'].values,
                                                                   fluence=qual_data_frame['fluence'].values,
                                                                   remaining_factor=qual_data_frame.values[:,-1],
                                                                   NIEL=data.proton_NIEL_SCREAM)

        elif args.particle is 'e':
            n_value = ddd_eq.fit_nValue_effectiveDDD(particle_energy=qual_data_frame['energy_mev'].values,
                                                     fluence=qual_data_frame['fluence'].values,
                                                     remaining_factor=qual_data_frame.values[:,-1],
                                                     NIEL=data.proton_NIEL_SCREAM,
                                                     energy_to_normalize=args.norm_e)

            ddd_vs_rf = ddd_eq.get_ddd_vs_remaining_factor(particle_energy=qual_data_frame['energy_mev'].values,
                                                                     fluence=qual_data_frame['fluence'].values,
                                                                     remaining_factor=qual_data_frame.values[:,-1],
                                                                     NIEL=data.proton_NIEL_SCREAM,
                                                                     energy_to_normalize=args.norm_e,
                                                                     n=n_value)

        if args.fit_type == 'single':
            fits = eq.singleDegradationEquation(ddd_vs_rf)
            coefficients = fits.fit(parameters=args.single_fit_param)

        elif args.fit_type == 'double':
            fits = eq.double_degradation_equation(ddd_vs_rf)
            coefficients = fits.fit(parameters=args.double_fit_param)

        else:
            fits = None


    # Calculate Remaining Factors
    if (args.remaining_factor is not None) and ((args.single_fit_param is not None) or (args.double_fit_param is not None)):
        if args.single_fit_param is not None:
            if len(args.single_fit_param) == 3:
                rf = eq.degradation_equation_rf_greater_than_one(fluence_ddd=args.remaining_factor, C=args.single_fit_param[0], phi_D_x=args.single_fit_param[1], A=args.single_fit_param[2])
            else:
                rf = eq.degradation_equation(fluence_ddd=args.remaining_factor, C=args.single_fit_param[0], phi_D_x=args.single_fit_param[1])
        elif args.double_fit_param is not None:
            if len(args.double_fit_param) == 5:
                rf = eq.double_degradation_equation_rf_greater_than_one(fluence_ddd=args.remaining_factor, C1=args.double_fit_param[0], phi_D_x1=args.double_fit_param[1], C2=args.double_fit_param[2], phi_D_x2=args.double_fit_param[3], A=args.double_fit_param[4])
            else:
                rf = eq.double_degradation_equation(fluence_ddd=args.remaining_factor, C1=args.double_fit_param[0], phi_D_x1=args.double_fit_param[1], C2=args.double_fit_param[2], phi_D_x2=args.double_fit_param[3])


    if args.output_json:
        if (ddd_vs_rf is not None) :
            ddd_vs_rf_dict = {'energy_mev':list(ddd_vs_rf[:,0]), 'ddd':list(ddd_vs_rf[:,1])}

            ddd_output_dict = {"particle":args.particle,
                               "param":qual_data_frame.columns[-1],
                               "qual_data": qual_data_frame.to_dict(orient='list'),
                               "ddd": ddd_vs_rf_dict}

            if coefficients is not None:
                fit_info = {}
                if args.fit_type == 'single':
                    fit_info = {'fit_type':'single', 'C':coefficients[0], 'phi':coefficients[1]}

                elif args.fit_type == 'double':
                    fit_info = {'fit_type':'double', 'C1':coefficients[0], 'phi1':coefficients[1], 'C2':coefficients[2], 'phi2':coefficients[3]}

                if fits is not None:
                    fit_info['r2'] = fits.r2

                ddd_output_dict['fit_parameters'] = fit_info

            # print(ddd_output_dict)
            os.chdir(args.output_loc)
            with open(args.output_json + '.json', 'w') as fp:
                json.dump(ddd_output_dict, fp, indent=4)

        if (args.remaining_factor is not None) and ((args.single_fit_param is not None) or (args.double_fit_param is not None)):
            rf_output = {'rf':rf}
            os.chdir(args.output_loc)
            # print('saving')
            with open('remaining_factor_ddd.json', 'w') as fp:
                json.dump(rf_output, fp, indent=4)

if __name__ == '__main__':
    print('This is your DDD calculator')
    main()
