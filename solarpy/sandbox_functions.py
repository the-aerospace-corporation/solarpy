import numpy as np
import scipy as sp
from scipy import interpolate
import solarpy.data_standards as sy_data
import solarpy.ddd.ddd_degradation_functions as sy_ddd_deq
import solarpy.eqflux as sy_eqflux

def adjust_niel(rdc_curve, niel_curve=None, normalizer=None, start_energy=None):
    if niel_curve is None:
        niel_curve = np.copy(sy_data.proton_NIEL_SCREAM)
    if normalizer is None:
        normalizer = 10
    rdc_curve = rdc_curve[rdc_curve[:,0]>=np.min(niel_curve[:,0])]
    rdc_curve = rdc_curve[rdc_curve[:,0]<=np.max(niel_curve[:,0])]
    rdc_curve = renorm_rdc(rdc_curve, energy_to_normalize=normalizer)
    ##### adjusting niel
    # interpolaate niel
    normalized_NIEL = niel_curve
    niel_normalizer = sy_ddd_deq.lookupNIELValue(normalizer, sy_data.proton_NIEL_SCREAM)
    normalized_NIEL[:,1] = niel_curve[:,1] / niel_normalizer
    intp_NIEL = sp.interpolate.interp1d(normalized_NIEL[:,0], normalized_NIEL[:,1])
    niel_intp_for_rdc = np.copy(rdc_curve)
    niel_intp_for_rdc[:, 1] = intp_NIEL(rdc_curve[:,0])

    rdc_below_start = rdc_curve[rdc_curve[:, 0] < start_energy]
    niel_below_Start = niel_intp_for_rdc[niel_intp_for_rdc[:, 0] < start_energy]
    # print(rdc_below_start[-1,0])
    # print(proton_gaas_NIEL_below_Start[-1,0])
    correction = rdc_below_start[:,1]/niel_below_Start[:,1]
    niel_below_Start[:,1] = (niel_below_Start[:,1]*correction)*niel_normalizer
    cut_niel = np.copy(sy_data.proton_NIEL_SCREAM[sy_data.proton_NIEL_SCREAM[:,0] >= start_energy])
    new_niel = np.vstack((niel_below_Start, cut_niel))
    return new_niel

def renorm_rdc(rdc_curve, energy_to_normalize):
    new_rdc = np.copy(rdc_curve)
    intp_func = sp.interpolate.interp1d(rdc_curve[:,0], rdc_curve[:,1])
    normalizer = intp_func(energy_to_normalize)
    new_rdc[:,1] = new_rdc[:,1]/normalizer
    return new_rdc

def fill_rdc(rdc_curve, filler_rdc_curve, minimum_particle_energy=1e-5, maximum_particle_energy=1e2):
    rdc_filler = np.vstack((filler_rdc_curve[filler_rdc_curve[:, 0] < np.min(rdc_curve[:, 0])], rdc_curve))
    rdc_filler = np.vstack((rdc_filler, filler_rdc_curve[filler_rdc_curve[:, 0] > np.max(rdc_curve[:, 0])]))
    rdc_filler = sy_eqflux.extrapolate_RDC_loglog(rdc_filler, minimum_particle_energy=minimum_particle_energy,maximum_particle_energy=maximum_particle_energy)
    return rdc_filler




