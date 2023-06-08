import argparse
import pandas as pd
import json
import os
from solarpy import data
from solarpy import eqflux
from solarpy import file_imports as sy_imports
import solarpy.degradation_equations as eq


electron_rdc_vars = ['electron qual data', 'electrical parameter', 'critical remaining factor',
                     'energy to normalize rdc', 'shield thickness in cm', 'density of shield', 'shield range table']
proton_rdc_vars = ['proton qual data', 'electrical parameter', 'critical remaining factor', 'energy to normalize rdc',
                   'shield thickness in cm', 'solar cell thickness' 'density of shield', 'shield range table',
                   'solar cell range table', 'solar cell displacement table']
params = ['pmax', 'voc', 'vmax', 'isc', 'imax', 'fillfactor', 'efficiency']


def main():
    parser = argparse.ArgumentParser(prog='RDC Calculator')
    parser.add_argument('-p', '--particle',
                        help='Calculate RDC for particle type: protons (p), electrons (e), both (b)',
                        choices=['p', 'e', 'b'])
    parser.add_argument('-param', help='Electrical Parameter', choices=params)
    parser.add_argument('-f', '--fit_type',
                        help='Select fit type: Single Degradation Equation (s), Double Degradation Equation(d)',
                        choices=['single', 'double'], default='single')
    parser.add_argument('-sfp', '--single_fit_param', help='Starting guess to fit C and phi or to calculate remaining factor', default=None, nargs=2, type=float, metavar='C Phi')
    parser.add_argument('-dfp', '--double_fit_param', help='Starting guess to fit C1, phi1, C2, phi2 or to calculate remaining factor', default=None, nargs=4, type=float, metavar='C1 phi1 C2 phi2')
    parser.add_argument('-fr', '--fit_rad_data', help='2d data of fluence and remaining factor', metavar='')
    parser.add_argument('-rf', '--remaining_factor', help='Get remaining factor by providing fluence', default = None, metavar='', type=float)
    parser.add_argument('-cf', '--conversion_factor', help='Get conversion factor from protons to electrons or vice versa.', metavar='rf | C_p phi_p | C_e phi_e, or rf | C_p1 phi_p1 C_p2 phi_p2 | C_e1 phi_e1 C_e2 phi_e2', nargs='*', type = float) # calculate conversion factors
    parser.add_argument('-i_rdc_u', '--input_unidirectional_rdc', help='Input unidirectional RDC curve to convert to shielded RDC')
    parser.add_argument('-r', '--rad', help='Radiation Data in either xlsx or txt', metavar='')
    parser.add_argument('-crf_p', help='critical remaining factor for protons', type=float, default=0.8, metavar='')
    parser.add_argument('-crf_e', help='critical remaining factor for electrons', type=float, default=0.8, metavar='')
    parser.add_argument('-norm_e', help='Electron energy to normalize electron radiation data', type=float, default=1,
                        metavar='')
    parser.add_argument('-norm_p', help='Proton energy to normalize proton radiation data', type=float, default=10,
                        metavar='')
    parser.add_argument('-st', help='Shield thickness in cm', type=float, default=0, metavar='')
    parser.add_argument('-sd', help='Shield density in g/cm2', type=float, default=2.554, metavar='')
    parser.add_argument('-sct', help='Solar cell thickness in cm', type=float, default=1e-2, metavar='')
    parser.add_argument('-srt', help='Shield range table xls', metavar='')
    parser.add_argument('-scrt', help='Solar cell range table xls', metavar='')
    parser.add_argument('-scdt', help='Solar cell displacement table xls', metavar='')
    parser.add_argument('-rdc', help='Return RDC curves: Unidirectional (u), Unidirectional Interpolated (ui), Omnidirectional Shielded (os) and All (a)',
                        choices=['u', 'ui', 'os', 'a'], default='a')
    parser.add_argument('-od', '--output_loc', help='Where to put output files', metavar='')
    parser.add_argument('-o', '--output_xlsx', help='Output Table. Give filename', metavar='')
    parser.add_argument('-oj', '--output_json', help='Output json. Give filename', metavar='')
    parser.add_argument('-oa', '--output_all', help='Output json. Give filename', metavar='')



    args = parser.parse_args()

    #initialize data lists
    coefficients = []
    fits = []
    rf = []
    cf = []

    if args.scrt:
        # TODO: Advanced feature to add custom range table
        print('Feature not enabled. Uses default GaAs range table from The GaAs Radiation Handbook')

    else:
        solar_cell_proton_range_table = data.GaAs_proton_ranges

    if args.scdt:
        # TODO: Advanced feature to add custom range table
        print('Feature not enabled. Uses default GaAs displacement table that was generated using TRIM/SRIM 2013')
    else:
        solar_cell_displacement_table = data.GaAs_total_displacements

    if args.srt:
        # TODO: Advanced feature to add custom range table
        print('Feature not enabled. Uses CMG proton and electron shield tables from The GaAs Radiation Handbook')

    else:
        shield_proton_range_table = data.CMG_proton_ranges
        shield_electron_range_table = data.CMG_electron_ranges

    if args.rad:
        qual_data_frame = sy_imports.get_dataframe_from_input_file(args.rad, type='qual')
        # print(qual_data_frame)
        if args.fit_type == 'single':
            fit_parameters = args.single_fit_param

        elif args.fit_type == 'double':
            fit_parameters = args.double_fit_param

        if args.particle == 'e':
                rdc_object = eqflux.electron_rdc_aero(qual_data=qual_data_frame.values,
                                                      critical_remaining_factor=args.crf_e,
                                                      energy_to_normalize_rdc=args.norm_e,
                                                      fit_type=args.fit_type,
                                                      fit_parameters=fit_parameters,
                                                      density_of_shield=args.sd,
                                                      shield_thickness_cm=args.st,
                                                      shield_range_table=shield_electron_range_table)

                fits = eqflux.reference_particle_fit(qual_data=qual_data_frame.values,
                                                     reference_energy=args.norm_e,
                                                     fit_type=args.fit_type,
                                                     fit_parameters=fit_parameters,
                                                     zero_fit=None
                                                     )

                coefficients = fits.coefficients

        elif args.particle == 'p':
            rdc_object = eqflux.proton_rdc_aero(qual_data=qual_data_frame.values, critical_remaining_factor=args.crf_p,
                                                energy_to_normalize_rdc=args.norm_p,
                                                fit_type=args.fit_type,
                                                fit_parameters=fit_parameters,
                                                density_of_shield=args.sd,
                                                shield_thickness_cm=args.st,
                                                shield_range_table=shield_proton_range_table,
                                                solar_cell_range_table=solar_cell_proton_range_table,
                                                solar_cell_displacement_table=solar_cell_displacement_table,
                                                solar_cell_thickness_cm=args.sct)

            fits = eqflux.reference_particle_fit(qual_data=qual_data_frame.values,
                                                 reference_energy=args.norm_p,
                                                 fit_type=args.fit_type,
                                                 fit_parameters=fit_parameters,
                                                 zero_fit=None)

            coefficients = fits.coefficients


        else:
            rdc_object = None

    # Takes user supplied unidirection RDC and calculates omnidirectional and/or shielded RDC
    if args.input_unidirectional_rdc:
        rdc_import = sy_imports.get_dataframe_from_input_file(args.input_rdc_unidirectional)

        if args.particle == 'e':
            rdc_shielded = eqflux.get_omnidirectional_electron_RDC_with_shielding(energy_vs_rdc=rdc_import.values,
                                                                                  mass_thickness_of_shield=args.sd * args.st,
                                                                                  shield_range_table=shield_electron_range_table
                                                                                  )

        elif args.particle == 'p':
            rdc_shielded = eqflux.get_omnidirectional_and_shielded_proton_rdc(energy_vs_rdc=rdc_import.proton_rdc,
                                                                            solar_cell_range_table=solar_cell_proton_range_table,
                                                                            shield_range_table=shield_proton_range_table,
                                                                            solar_cell_displacement_table=solar_cell_displacement_table,
                                                                            solar_cell_thickness_cm=args.sct,
                                                                            shield_thickness_cm=args.st
                                                                            )
        else:
            rdc_shielded = None


    if args.fit_rad_data:
        rad_data = sy_imports.get_dataframe_from_input_file(args.fit_rad_data)

        if args.fit_type == 'single':
            fit = eq.singleDegradationEquation(fluence_vs_remaining_factor=rad_data.values)
            coefficients = fit.fit(parameters=args.single_fit_param)

        elif args.fit_type == 'double':
            fit = eq.double_degradation_equation(fluence_vs_remaining_factor=rad_data.values)
            coefficients = fit.fit(parameters=args.double_fit_param)

    # Calculate Remaining Factors
    if (args.remaining_factor is not None) and ((args.single_fit_param is not None) or (args.double_fit_param is not None)):
        if args.single_fit_param is not None:
            rf = eq.degradation_equation(fluence_ddd=args.remaining_factor, C=args.single_fit_param[0], phi_D_x=args.single_fit_param[1])
        elif args.double_fit_param is not None:
            rf = eq.double_degradation_equation(fluence_ddd=args.remaining_factor, C1=args.double_fit_param[0], phi_D_x1=args.double_fit_param[1], C2=args.double_fit_param[2], phi_D_x2=args.double_fit_param[3])

    # Calculate Proton to Electron or Electron to Proton Conversion
    if args.conversion_factor:
        if len(args.conversion_factor <= 5):
            rf_p = eq.get_fluence(remaining_factor=args.conversion_factor[0],
                                  C=args.conversion_factor[1],
                                  phi_D_x=args.conversion_factor[2])

            rf_e = eq.get_fluence(remaining_factor=args.conversion_factor[0],
                                  C=args.conversion_factor[3],
                                  phi_D_x=args.conversion_factor[4])

            cf = rf_e/rf_p

        elif len(args.conversion_factor > 5):
            rf_p = eq.get_fluence_ddd_double_degradation(remaining_factor=args.conversion_factor[0],
                                                         C1=args.conversion_factor[1],
                                                         phi_D_x1=args.conversion_factor[2],
                                                         C2=args.conversion_factor[3],
                                                         phi_D_x2=args.conversion_factor[4])

            rf_e = eq.get_fluence_ddd_double_degradation(remaining_factor=args.conversion_factor[0],
                                                         C1=args.conversion_factor[5],
                                                         phi_D_x1=args.conversion_factor[6],
                                                         C2=args.conversion_factor[7],
                                                         phi_D_x2=args.conversion_factor[8])

            cf = rf_e / rf_p

    if args.output_json:
        if args.rad:
            if args.rdc == 'u':
                rdc = {"unidirectional": {"energy_mev": list(rdc_object.energy_vs_rdc[:, 0]),"rdc": list(rdc_object.energy_vs_rdc[:, 1])}}

            elif args.rdc == 'ui':
                rdc = {"unidirectional_extrapolated": {"energy_mev": list(rdc_object.energy_vs_rdc_extrapolated[:, 0]),"rdc": list(rdc_object.energy_vs_rdc_extrapolated[:, 1])}}

            elif args.rdc == 'os':
                rdc = {"omnidirectional_shield": {"energy_mev": list(rdc_object.omnidirectional_shielded_rdc[:, 0]),"rdc": list(rdc_object.omnidirectional_shielded_rdc[:, 1])}}

            else:
                rdc= {"unidirectional": {"energy_mev": list(rdc_object.energy_vs_rdc[:,0]), "rdc": list(rdc_object.energy_vs_rdc[:,1])},
                      "unidirectional_extrapolated": {"energy_mev": list(rdc_object.energy_vs_rdc_extrapolated[:,0]), "rdc": list(rdc_object.energy_vs_rdc_extrapolated[:,1])},
                      "omnidirectional_shield": {"energy_mev": list(rdc_object.omnidirectional_shielded_rdc[:,0]), "rdc": list(rdc_object.omnidirectional_shielded_rdc[:,1])}
                      }

            if coefficients is not None:
                fit_info = {}
                if args.fit_type == 'single':
                    fit_info = {'fit_type':'single', 'C':coefficients[0], 'phi':coefficients[1]}

                elif args.fit_type == 'double':
                    fit_info = {'fit_type':'double', 'C1':coefficients[0], 'phi1':coefficients[1], 'C2':coefficients[2], 'phi2':coefficients[3]}

                if fits is not None:
                    fit_info['r2'] = fits.r2

            output_dict = {"particle":args.particle, "param":qual_data_frame.columns[-1], "qual_data": qual_data_frame.to_dict(orient='list'), "reference_energy":rdc_object.energy_to_normalize_rdc,"fit_parameters":fit_info ,"rdc": rdc}

        elif args.input_unidirectional_rdc:
            if args.rdc == 'ui':
                rdc = {"unidirectional_extrapolated": {"energy_mev": list(rdc_import.values[:, 0]),"rdc": list(rdc_import.values[:, 1])}}

            elif args.rdc == 'os':
                rdc = {"omnidirectional_shield": {"energy_mev": list(rdc_shielded[:,0]), "rdc": list(rdc_shielded[:,1])}}

            else:
                rdc = {"unidirectional_extrapolated": {"energy_mev": list(rdc_import.values[:, 0]),"rdc": list(rdc_import.values[:, 1])},
                       "omnidirectional_shield": {"energy_mev": list(rdc_shielded[:,0]), "rdc": list(rdc_shielded[:,1])}
                       }
            output_dict = {"particle":args.particle, "param":'fill in', "rdc": rdc}

        os.chdir(args.output_loc)
        with open(args.output_json+'.json', 'w') as fp:
            json.dump(output_dict,fp, indent=4)


            # # os.chdir('C:\Users\dw28596\PycharmProjects\solarpy CLI')
            # output_file = OutputCellTechRDCdata(electron_rdc=electron_rdc,
            #                                     proton_rdc=proton_rdc,
            #                                     type_of_rdc=args.rdc,
            #                                     params='All',
            #                                     )
            # output_file.output(args.output_xlsx)

        if (args.remaining_factor is not None) and ((args.single_fit_param is not None) or (args.double_fit_param is not None)):
            rf_output = {'rf':rf}
            os.chdir(args.output_loc)
            # print('saving')
            with open('remaining_factor_eqflux.json', 'w') as fp:
                json.dump(rf_output, fp, indent=4)

        if args.conversion_factor:
            cf_output = {'cf':cf}
            os.chdir(args.output_loc)
            # print('saving')
            with open('conversion_factor.json', 'w') as fp:
                json.dump(cf_output, fp, indent=4)



    if args.output_xlsx:
        print('output xlsx')
        filename = args.output_xlsx + '.xlsx'
        os.chdir(args.output_loc)
        if args.rad:
            with pd.ExcelWriter(filename) as writer:
                qual_data_frame.to_excel(writer, sheet_name='qual_data', index=False)
                if args.rdc == 'u':
                    pd.DataFrame(data=rdc_object.energy_vs_rdc, columns=['energy_mev', 'rdc']).to_excel(writer, sheet_name='unidirectional', index=False)

                elif args.rdc == 'ui':
                    pd.DataFrame(data=rdc_object.energy_vs_rdc_extrapolated, columns=['energy_mev', 'rdc']).to_excel(writer, sheet_name='unidirectional_extrapolated', index=False)

                elif args.rdc == 'os':
                    pd.DataFrame(data=rdc_object.omnidirectional_shielded_rdc, columns=['energy_mev', str(args.st/0.00254)+' mil']).to_excel(writer, sheet_name='omnidirectional_shielded', index=False)

                else:
                    pd.DataFrame(data=rdc_object.energy_vs_rdc, columns=['energy_mev', 'rdc']).to_excel(writer, sheet_name='unidirectional', index=False)
                    pd.DataFrame(data=rdc_object.energy_vs_rdc_extrapolated, columns=['energy_mev', 'rdc']).to_excel(writer, sheet_name='unidirectional_extrapolated', index=False)
                    pd.DataFrame(data=rdc_object.omnidirectional_shielded_rdc, columns=['energy_mev', str(args.st/0.00254)+' mil']).to_excel(writer, sheet_name='omnidirectional_shielded', index=False)

        elif args.input_unidirectional_rdc:
            with pd.ExcelWriter(filename) as writer:
                if args.rdc == 'ui':
                    pd.DataFrame(data=rdc_import.values, columns=['energy_mev', 'rdc']).to_excel(writer, sheet_name='unidirectional_extrapolated',index=False)

                elif args.rdc == 'os':
                    pd.DataFrame(data=rdc_shielded, columns=['energy_mev', str(args.st / 0.00254) + ' mil']).to_excel(writer, sheet_name='omnidirectional_shielded', index=False)

                else:
                    pd.DataFrame(data=rdc_import.values, columns=['energy_mev', 'rdc']).to_excel(writer, sheet_name='unidirectional_extrapolated',index=False)
                    pd.DataFrame(data=rdc_shielded, columns=['energy_mev', str(args.st / 0.00254) + ' mil']).to_excel(writer, sheet_name='omnidirectional_shielded', index=False)

    #         # os.chdir(args.ol)
    #         output_file = OutputCellTechRDCdata(electron_rdc=rdc_object,
    #                                             proton_rdc=rdc_object,
    #                                             type_of_rdc=args.rdc,
    #                                             params='All',
    #                                             )
    #         output_file.output_json(args.output_json)

    # if args.output_all:
    #     print('output all')
    #     if args.rad:
    #         os.chdir(args.output_loc)
    #         output_file = OutputCellTechRDCdata(electron_rdc=rdc_object,
    #                                             proton_rdc=rdc_object,
    #                                             type_of_rdc=args.rdc,
    #                                             params='All',
    #                                             )
    #         output_file.output_all(args.output_all)


def rdc_arg_parser(arg, rdc_param):
    if arg == 'u':
        print(rdc_param.get_energy_vs_rdc())

    if arg =='ui':
        print(rdc_param.get_energy_vs_rdc_extrapolated())

    if arg == 'os':
        print(rdc_param.get_omnidirectional_shielded_rdc())

    if arg == 'a':
        print(rdc_param.get_energy_vs_rdc())
        print(rdc_param.get_energy_vs_rdc_extrapolated())
        print(rdc_param.get_omnidirectional_shielded_rdc())


if __name__ == '__main__':
    print('This is your RDC calculator')
    main()

