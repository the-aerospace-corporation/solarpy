import pandas as pd
import numpy as np


def get_dataframe_from_input_file(file_name, type=None):
    """
    Gets data from file regardless of file type and outputs it as a dataframe

    Args:
        file_name: a string that contains the full file path
        type: '.txt', '.json', '.xlsx'

    Returns:

    """
    if file_name.endswith('.txt'):
        data_frame = pd.read_csv(file_name, sep='\t')

    elif file_name.endswith('.json'):
        data_frame = pd.read_json(file_name)

    elif file_name.endswith('.xlsx'):
        data_frame = pd.ExcelFile(file_name)

    else:
        data_frame = None
        print('No SOUP FOR YOU!')

    if type == 'qual':
        data_frame = parse_and_organize_qual_dataframe(data_frame)

    if type == 'env':
        data_frame = parse_and_oragnize_env_dataframe(data_frame)

    return data_frame

def parse_and_organize_qual_dataframe(data_frame=pd.DataFrame()):
    new_dateframe = pd.DataFrame(columns=['energy_mev', 'fluence'])
    for column in data_frame.columns:
        if column.lower() in ['energy (mev)', 'energy mev', 'energy', 'energy_mev']:
            new_dateframe.energy_mev = data_frame[column]

        elif column.lower() in ['fluence', 'fluence p/cm2', 'fluence (p/cm2)', 'fluence(e/cm2)', 'fluence(p/cm2']:
            new_dateframe.fluence = data_frame[column]

        elif column.lower() in ['pmax', 'p max', 'p_max', 'npmp', 'pmp', 'p_mp', 'p mp']:
            new_dateframe['pmax'] = data_frame[column]

        elif column.lower() in ['efficiency', 'eff']:
            new_dateframe['efficiency'] = data_frame[column]

        elif column.lower() in ['isc', 'i sc', 'i_sc','nisc']:
            new_dateframe['isc'] = data_frame[column]

        elif column.lower() in ['voc', 'v oc', 'v_oc', 'nvoc']:
            new_dateframe['voc'] = data_frame[column]

        elif column.lower() in ['imax', 'i_max', 'i max', 'imp', 'i mp', 'i_mp','nimp']:
            new_dateframe['imax'] = data_frame[column]

        elif column.lower() in ['vmax', 'v_max', 'v max', 'vmp', 'v mp', 'v_mp','nvmp']:
            new_dateframe['vmax'] = data_frame[column]

        elif column.lower() in ['isc', 'i sc', 'i_isc','nisc']:
            new_dateframe['isc'] = data_frame[column]

    return new_dateframe

def parse_and_oragnize_env_dataframe(data_frame=pd.DataFrame()):
    env = environment()
    env.proton = np.vstack((data_frame['proton']['energy_mev'], data_frame['proton']['fluence'])).T
    env.electron = np.vstack((data_frame['electron']['energy_mev'], data_frame['electron']['fluence']))
    return env

class environment():
    def __init__(self):
        self.proton = []
        self.electron = []