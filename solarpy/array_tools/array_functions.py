import numpy as np
import scipy as sp
from scipy import interpolate

def current_loss_factors_linear(eol_fluence_time, coverglass_darkening=1, contamination=1, loss_from_random_failures=1):
    total_linear_loss_factor = coverglass_darkening*contamination*loss_from_random_failures
    slope_of_loss_factor = (1 - total_linear_loss_factor)/eol_fluence_time
    return slope_of_loss_factor