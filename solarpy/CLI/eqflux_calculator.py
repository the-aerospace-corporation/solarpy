import argparse
from solarpy import eqflux
from solarpy.file_imports.ImportRDCdata import ImportRDCdata
from solarpy.file_imports.ImportEnvironmentSpectra import ImportEnvironmentSpectra

params = ['Pmax', 'Voc', 'Vmax', 'Isc', 'Imax', 'FillFactor', 'Efficiency']


def main():
    parser = argparse.ArgumentParser(prog='EQFLUX Calculator')
    parser.add_argument('-rdc', help='RDC xlsx File')
    parser.add_argument('-env', help='Environment Spectra')
    parser.add_argument('-y', help='Number of years to multiply by', type=float, default=1, metavar='')

    args = parser.parse_args()
    print(args)

    if args.env:
        env = ImportEnvironmentSpectra(args.env)

    if args.input_rdc: # make this a class
        rdc = ImportRDCdata(args.input_rdc)
        for param in params:
            if getattr(rdc.electron_rdc_interpolated, param):
                proton_rdc_param = getattr(rdc.proton_rdc_interpolated, param)
                ten_MeV_protons_differential = eqflux.get_relative_MeV_fluence_trapezoidal_integration(proton_rdc_param, env.protons_differential)
            else:
                ten_MeV_protons_differential = 0

            if getattr(rdc.proton_rdc_interpolated, param):
                electron_rdc_parm = getattr(rdc.electron_rdc_interpolated, param)
                one_MeV_electrons_differential = eqflux.get_relative_MeV_fluence_trapezoidal_integration(electron_rdc_parm, env.electrons_differential)
            else:
                one_MeV_electrons_differential = 0

            total_one_MeV_electrons_differential = one_MeV_electrons_differential + (ten_MeV_protons_differential * rdc.proton_to_electron_conv) ### need to work on this



if __name__ == '__main__':
    print('This is your EQFLUX calculator')
    main()
