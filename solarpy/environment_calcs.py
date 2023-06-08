import numpy as np
import pandas as pd
import solarpy.eqflux as eqf

def get_environment_over_time(env, start_time, period, unit):
    """
    Works with the flux output from Spenvis

    Args:
        env: 2d data with energy, flux (p/s), fluence
        start_time: datetime object of start
        period: how many units you want
        unit: units you want the time to be in m=minutes, h=hours, d=days, and s=seconds

    Returns:

    """
    repeats_time_multiplied = get_environment_per_time(env, start_time, period, unit)
    cumulative_fluence_dataframe = repeats_time_multiplied.cumsum()
    return cumulative_fluence_dataframe

def get_environment_per_time(env, start_time, period, unit):
    """
    Works with the flux output from Spenvis

    Args:
        env:
        start_time:
        period:
        unit:

    Returns:

    """
    if unit == 'm':
        time_multiplier = 60
        freq = 'min'
    elif unit == 'h':
        time_multiplier = 60*60
        freq = 'H'
    elif unit == 'd':
        time_multiplier = 24*60*60
        freq = 'D'
    else:
        time_multiplier = 1
        freq = 'S'

    repeats = np.repeat(np.array([env[:,1]]), period, axis=0)
    t_o = np.zeros((np.shape(env[:,1])))
    repeats = np.vstack((t_o,repeats))
    repeats_time_multiplied = repeats*time_multiplier # the env spectra is usually in particles/s, hence the multiplier when you want data over anything other seconds time steps
    fluence_period = pd.DataFrame(data=repeats_time_multiplied, columns=env[:,0], index=pd.date_range(start=start_time, periods=period+1, freq=freq))
    return fluence_period

def get_one_MeV_over_time(datetime, rdc, list_of_environments, method='aero', particle='e'):
    one_MeV_fluence_over_time = []
    for env in list_of_environments:
        if method == 'EQFLUX':
            if particle == 'e':
                one_MeV_fluence_over_time.append(eqf.compute_equivalent_fluence_for_electrons(rdc, env))
            elif particle == 'p':
                one_MeV_fluence_over_time.append(eqf.compute_equivalent_fluence_for_protons(rdc, env))
        else:
            one_MeV_fluence_over_time.append(eqf.get_relative_MeV_fluence_trapezoidal_integration(rdc, env))

    dataframe = pd.DataFrame(one_MeV_fluence_over_time, columns=['one_mev_fluence'])
    dataframe = dataframe.set_index(datetime)
    return dataframe

def one_mev_vs_time_list(dataframe):
    energy_mev = dataframe.columns
    mission_env_array = dataframe.to_numpy()
    mission_env_list = []
    for row in mission_env_array:
        mission_env_list.append(np.vstack((energy_mev, row)).astype(float).T)
    return mission_env_list

def get_mission_env_cum_sum_over_time(transfer_orbit_rad=None, orbit_rad=None, solar_protons=None, transfer_start_date=None, time_in_transfer=None, orbit_start_time=None, time_in_orbit=None):
    """
    combines radiation data from transfer orbit, orbit, and solar electric protons

    Args:
        transfer_orbit_rad:
        orbit_rad:
        transfer_start_date:
        time_in_transfer:
        orbit_start_time:
        time_in_orbit:
        solar_protons:

    Returns:

    """
    if transfer_orbit_rad is not None:
        gto_per_time = get_environment_per_time(transfer_orbit_rad, start_time=transfer_start_date, period=time_in_transfer, unit='d')
        geo_per_time = get_environment_per_time(orbit_rad, start_time=orbit_start_time, period=time_in_orbit, unit='d')
        mission_env_time = pd.concat([gto_per_time, geo_per_time], sort=True)

        if solar_protons is not None:
            sep_per_time = get_environment_per_time(solar_protons, start_time=transfer_start_date, period=time_in_orbit + 9, unit='d')
            mission_env = mission_env_time.add(sep_per_time, fill_value=0).cumsum()

        else:
            mission_env = pd.concat([gto_per_time, geo_per_time], sort=True).cumsum()

    else:
        geo_per_time = get_environment_per_time(orbit_rad, start_time=orbit_start_time, period=time_in_orbit, unit='d')
        mission_env_time = geo_per_time

        if solar_protons is not None:
            sep_per_time = get_environment_per_time(solar_protons, start_time=orbit_start_time, period=time_in_orbit, unit='d')
            mission_env = mission_env_time.add(sep_per_time, fill_value=0).cumsum()

        else:
            mission_env = mission_env_time.cumsum()

    return mission_env

