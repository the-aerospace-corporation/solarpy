import scipy as sp
import scipy.optimize
from scipy import interpolate
from scipy import integrate
import numpy as np
from solarpy import degradation_equations as eq
from solarpy.analyze_goodness_of_fit import analyze_goodness_of_fit
from solarpy import data_standards as data
import matplotlib.pyplot as plt

class DDDqualData:
    def __init__(self, qual_data, particle_type, NIEL=None, energy_to_normalize=None, fit_type=None, fit_parameters=None):
        self.qual_data = qual_data
        self.particle_type = particle_type
        self.particle_energy = qual_data[:,0]
        self.fluence = qual_data[:,1]
        self.remaining_factor = qual_data[:,2]
        self.particle_type

        if NIEL is None:
            if self.particle_type == 'e':
                self.NIEL = data.electron_NIEL_SCREAM
            elif self.particle_type == 'p':
                self.NIEL = data.proton_NIEL_SCREAM
        else:
            self.NIEL = NIEL

        if energy_to_normalize == None:
            self.energy_to_normalize = 1
        else:
            self.energy_to_normalize = energy_to_normalize

        if fit_type == None:
            self.fit_type = 'single'
        else:
            self.fit_type = fit_type

        if fit_parameters  == None:
            if self.fit_type == 'single':
                self.fit_parameters = None
            else:
                self.fit_parameters = fit_parameters
        else:
            if self.fit_type == 'double':
                self.fit_parameters = None
            else:
                self.fit_parameters = fit_parameters

        self.update_ddd()

    def update_ddd(self):
        if self.particle_type == 'e':
            self.n_value = fit_nValue_effectiveDDD(self.particle_energy, self.fluence, self.remaining_factor, self.NIEL, self.energy_to_normalize, self.fit_type)
        else:
            self.n_value = 1

        self.ddd_vs_rf = get_ddd_vs_remaining_factor(self.particle_energy, self.fluence, self.remaining_factor, self.NIEL, self.energy_to_normalize, self.n_value)

        if self.fit_type == 'single':
            self.fit = eq.singleDegradationEquation(self.ddd_vs_rf)
            self.fit.fit(parameters=self.fit_parameters)
        elif self.fit_type == 'double':
            self.fit = eq.doubleDegradationEquation(self.ddd_vs_rf)
            self.fit.fit(parameters=self.fit_parameters)

        self.coefficients = self.fit.coefficients
        self.r2 = self.fit.r2

    def plot_ddd_fits(self):
        by_energy = _groupByEnergy(self.particle_energy, self.qual_data)
        for energy in by_energy:
            ddd_vs_rf = get_ddd_vs_remaining_factor(energy[:,0], energy[:,1], energy[:,2], self.NIEL, self.energy_to_normalize, self.n_value)
            plt.plot(ddd_vs_rf[:,0], ddd_vs_rf[:,1], 'o', label=energy[0,0])
            plt.legend()

        x_ddd = np.logspace(np.log10(self.ddd_vs_rf[0,0]), np.log10(self.ddd_vs_rf[-1,0]), 50)
        if self.fit_type == 'single':
            rf_fitted = eq.degradation_equation(x_ddd, self.coefficients[0], self.coefficients[1])

        elif self.fit_type == 'double':
            rf_fitted = eq.double_degradation_equation(x_ddd, self.coefficients[0], self.coefficients[1], self.coefficients[2], self.coefficients[3])

        plt.semilogx(x_ddd, rf_fitted, '-', label='fit')
# def degradation_equation(D_d, C, D_x):
#     """
#     Degradation equation for determining remaining factor as a function of displacement damage dose.  The equation is exactly the same as the equation for EQFLUX
#
#     Args:
#         D_d (ndarray) : 1d numpy array of displacement damage dose
#         C (float): Constant that is determined through fitting
#         D_x (float) : displacement damage dose constant that is determined through fitting
#
#     Returns:
#         1d array of calculated remaining factors
#
#     """
#     D_x = np.abs(D_x)  # prevents negative numbers into log
#     remainingFactor = 1 - C * np.log10(1 + (D_d / D_x))
#     return remainingFactor
#
#
# def double_degradation_equation(D_d, C1, D_x1, C2, D_x2):
#     """
#     Combines two of the single degradation equation in an effort to alleviate the assumption that there are two components degrading.  Can fit more complex solar cell degradation better, but may not always be the aright approach as as it is an assumption that there are two mechanisms behaving similarly
#
#     Args:
#         D_d (ndarray) : 1d numpy array of displacement damage dose
#         C1 (float) : Constant that is determined through fitting
#         D_x1 (float) : displacement damage dose constant that is determined through fitting
#         C2 (float) : Constant that is determined through fitting
#         D_x2 (float) : displacement damage dose constant that is determined through fitting
#
#     Returns:
#         1d array of calculated remaining factors
#
#     """
#     D_x1 = np.abs(D_x1)  # prevents negative numbers into log
#     D_x2 = np.abs(D_x2)
#     remainingFactor = 1 - (C1 * np.log10(1 + (D_d / D_x1))) + (C2 * np.log10(1 + (D_d / D_x2)))
#     return remainingFactor
#
#
# def degradation_equation_RFgreater_than_one(D_d, A, C, D_x):
#     """
#     Degradation equation for determining remaining factor as a function of displacement damage dose.  The 1 in the degradation equation is replaced with a fitting parameter (A) to account for values above one.  This should not be used because it assumes that the remaining factors can be greater than beginning of life. It is in included as a reference to previous literature.
#
#     Args:
#         D_d (ndarray) : 1d numpy array of displacement damage dose
#         A (float) : Const that is determined from fitting to account for remaining factors greater than 1
#         C (float) : Constant that is determined through fitting
#         D_x (float) : displacement damage dose constant that is determined through fitting
#
#     Returns:
#         1d array of calculated remaining factors
#
#     """
#     D_x = np.abs(D_x)  # prevents negative numbers into log
#     remainingFactor = A - C * np.log10(1 + (D_d / D_x))
#     return remainingFactor
#
#
def lookupNIELValue(particle_energy, NIEL):
    """
    Finds a NIEL for a given particle energy.  Basically interpolates the data to find the NIEL

    Args:
        particle_energy: particle energy you want the NIEL value for
        NIEL (ndarray) : the 2D numpy array of particle energy vs NIEL

    Returns:
        NIEL value at the chosen particle energy

    """
    intp_data = sp.interpolate.interp1d(NIEL[:, 0], NIEL[:, 1], kind='linear', fill_value = 'extrapolate')

    # if np.size(particle_energy) > 1:
    #     max_energy = np.max(particle_energy)
    #     min_energy = np.min(particle_energy)
    # else:
    #     max_energy = particle_energy
    #     min_energy = particle_energy
    #
    # if (max_energy > np.max(NIEL[:, 0])) or (min_energy < np.min(NIEL[:, 0])):
    #     foundNIEL = 0
    # else:
    #     foundNIEL = intp_data(particle_energy)
    foundNIEL = intp_data(particle_energy)
    return foundNIEL


def get_ddd(particle_energy, fluence, NIEL, energy_to_normalize=None, n_value=None):
    """
    Calculates the  displacement damage dose.  For protons the "n" value is 1 and the energy to normalize isn't necessary. When calculating effective displacement damage dose an empirically determined "n" value is used to collapse the electron displacement damage curves for different particle energies.  The n value can be derived by fitting displacement damage curves to for a n until they collapse.  N can be anywhere from 0 to 4 as determined in the literature

    Args:
        particle_energy: particle energy
        fluence: fluence
        NIEL: 2d numpy array of particle energies in column 0 and nonionizing energy loss (NIEL) in column 1
        energy_to_normalize: particle energy to normalize the effective displacement damage dose curve. This value is typically 1 for 1MeV electrons
        n_value: Empirically determined fitting parameter to collapse displacement damage dose curves that do not collapse. It is 1 for protons

    Returns:
        1d numpy array of effective displacement damage dose

    """
    if energy_to_normalize == None:
        normalizedNIEL = 1
    else:
        normalizedNIEL = lookupNIELValue(energy_to_normalize, NIEL)

    if n_value == None:
        n_value = 1

    NIEL_value = lookupNIELValue(particle_energy, NIEL)
    displacementDamageDose = fluence * NIEL_value * ((NIEL_value / normalizedNIEL) ** (n_value - 1))
    D_d = displacementDamageDose

    return np.array(D_d)


def get_ddd_vs_remaining_factor(particle_energy, fluence, remaining_factor, NIEL, energy_to_normalize=None, n=None):
    """
    Calculates the 2d numpy array of displacement damage dose vs energy.

    Args:
        particle_energy (ndarray or float): particle energy
        fluence (ndarray) : 1d numpy array of fluence
        remaining_factor (ndarray) : 1d array of remaining factor of interest
        NIEL (ndarray): 2d numpy array of particle energies in column 0 and nonionizing energy loss (NIEL) in column 1
        energy_to_normalize (float) : particle energy to normalize the effective displacement damage dose curve. This value is typically 1 for 1MeV electrons
        n (float): Empirically determined fitting parameter to collapse displacement damage dose curves that do not collapse. It is 1 for protons

    Returns:
        2d numpy array where column 0 is the displacement damage dose and column 1 is the remaining factor

    """
    D_e = get_ddd(particle_energy, fluence, NIEL, energy_to_normalize, n)
    DDDvsRemainingFactor = np.vstack((D_e, remaining_factor)).transpose()
    DDDvsRemainingFactor = DDDvsRemainingFactor[np.argsort(DDDvsRemainingFactor[:, 0])]
    return DDDvsRemainingFactor


def error_function_nValue_for_effectiveDDD(parameter, particle_energy, fluence, remaining_factor, NIEL,
                                           energy_to_normalize, fit_type):
    """
    Error function used to find the the n value for radiation qualification data that does not collapse when using displacement damage dose.  For protons this has empirically determined to be 1 and for electrons the n value can range from 0.3-4

    Args:
        parameter: n value guess to minimize
        particle_energy: particle energy
        fluence:  1d numpy array of fluence
        remaining_factor: 1d array of remaining factor of interest
        NIEL: 2d numpy array of particle energies in column 0 and nonionizing energy loss (NIEL) in column 1
        energy_to_normalize: particle energy to normalize the effective displacement damage dose curve. This value is typically 1 for 1MeV electrons
        fit_type: fit using single degradation equation, 'single', or double degradation equation, 'double'.  May add another option to fit using the A parameter

    Returns:
        Error to be minimized by using 1 - r2.  # TODO: will look into using chi_squared instead
    """
    DDDvsRF = get_ddd_vs_remaining_factor(particle_energy, fluence, remaining_factor, NIEL,
                                          energy_to_normalize, parameter)
    if fit_type == 'single':
        fit = eq.fit_degradation_equation(DDDvsRF, parameters=[2.5e-2, 2.0e8])
        newRF = eq.degradation_equation(DDDvsRF[:, 0], fit[0], fit[1])
    elif fit_type == 'double':
        fit = eq.fit_double_degradation_equation(DDDvsRF, parameters=[1e-2, 1e2, 1e2, 1e-2]) # [1e-2, 1e2, 1e3, 1e-2]
        newRF = eq.double_degradation_equation(DDDvsRF[:, 0], fit[0], fit[1], fit[2], fit[3])
    goodnessOfFit = analyze_goodness_of_fit(DDDvsRF[:, 1], newRF)
    error = 1 - goodnessOfFit.r2
    return error


# def error_function_degradation_equation(p, D_d, remainingFactor):
#     numberDDDpoints = len(D_d)
#     RF_new = (((np.abs(remainingFactor - degradation_equation(D_d, *p)) / remainingFactor) * 100) ** 2)
#     error = np.sqrt(RF_new.sum() / numberDDDpoints)
#     return error
#
#
# def error_function_double_degradation_equation(p, D_d, remainingFactor):
#     numberDDDpoints = len(D_d)
#     RF_new = (((np.abs(remainingFactor - double_degradation_equation(D_d, *p)) / remainingFactor) * 100) ** 2)
#     error = np.sqrt(RF_new.sum() / numberDDDpoints)
#     return error
#
#
# def error_function_degradation_equation_RFgreater_than_one(p, D_d, remainingFactor):
#     numberDDDpoints = len(D_d)
#     RF_new = (((np.abs(
#         remainingFactor - degradation_equation_RFgreater_than_one(D_d, *p)) / remainingFactor) * 100) ** 2)
#     error = np.sqrt(RF_new.sum() / numberDDDpoints)
#     return error
#
#
# def fitDegradationEquation(DDD_vs_remaining_factor, parameters=None):
#     """
#     Given the displacement damage dose curve and minimizing using the Nelder-Mead method to arrive at the best fit for C and D_x for the DDD degradation equation.  The data is fit by minimizing the error function to the root mean squared error (RMSE).  This is done becuase it gives lower chi squared values and r squared values closer to 1.  Also, you can be pretty far off on your starting parameters and still arrive at a good fit.  RMSE allows for better fits as there is usually error in the fluence and remaining factor and provides the best fit to expected values.
#
#     Args:
#         DDD_vs_remaining_factor (ndarray) : 2d numpy array where column 0 is the displacement damage dose and column 1 is the remaining factor
#         parameters (ndarray): 1d numpy array of starting values for the fit. The array should be 2 elements where the first element is the C parameter in the degradation equation and D_x is the second parameter. If no parameters are entered the fit defaults to using 2e-2 for C and 2e8 for D_x.  These parameters are found to fit most data without any modification
#
#     Returns:
#         1d numpy array of 2 elements where the first element is the C parameter and the second is D_x
#
#     """
#     if parameters == None:
#         parameters = [2e-2, 2e8]
#
#     if DDD_vs_remaining_factor[0, 0] != 0:
#         D_d = np.vstack(([0, 1], DDD_vs_remaining_factor))  # make sure we fit to 0 and 1
#     else:
#         D_d = DDD_vs_remaining_factor
#
#     fit = sp.optimize.minimize(error_function_degradation_equation, parameters, args=(D_d[:, 0], D_d[:, 1]),
#                                method='Nelder-Mead',
#                                options={'maxiter': 10000, 'maxfev': 10000, 'ftol': 1e-6, 'xtol': 1e-6})  # tol=1e-8) #
#     fit.x = np.abs(fit.x)  # deals with negative parameters for log
#     return fit.x
#
#
# def fitDoubleDegradationEquation(DDD_vs_remaining_factor, parameters=None):
#     """
#     Given the displacement damage dose curve and minimizing using the Nelder-Mead method to arrive at the best fit for C1, D_x1, C2, D_x2 for the DDD double degradation equation.  The data is fit by minimizing the error function to the root mean squared error (RMSE).  This is done becuase it gives lower chi squared values and r squared values closer to 1.  Also, you can be pretty far off on your starting parameters and still arrive at a good fit.  RMSE allows for better fits as there is usually error in the fluence and remaining factor and provides the best fit to expected values.
#
#     Args:
#         DDD_vs_remaining_factor (ndarray) : 2d numpy array where column 0 is the displacement damage dose and column 1 is the remaining factor
#         parameters (ndarray): 1d numpy array of starting values for the fit. The array should be 4 elements where the first element is the C1 parameter in the degradation equation, D_x1 is the second parameter, C2 is the third parameter, and D_x2 is the last paramenter.  If no parameters are entered the fit defaults to using 6e-2, 6e8, 2e-4, 4e-10 for C1, D_x1, C2, and D_x2 respectively.  These parameters are found to fit most data without any modification
#     Returns:
#         1d numpy array of 4 elements where the first element is the C1 parameter in the degradation equation, D_x1 is the second parameter, C2 is the third parameter, and D_x2 is the last paramenter.
#
#     """
#     if parameters == None:
#         parameters = [6e-2, 6e8, 2e-4, 4e-10]
#
#     if DDD_vs_remaining_factor[0, 0] != 0:
#         D_d = np.vstack(([0, 1], DDD_vs_remaining_factor))  # make sure we fit to 0 and 1
#     else:
#         D_d = DDD_vs_remaining_factor
#
#     fit = sp.optimize.minimize(error_function_double_degradation_equation, parameters, args=(D_d[:, 0], D_d[:, 1]),
#                                method='Nelder-Mead',
#                                options={'maxiter': 10000, 'maxfev': 10000, 'ftol': 1e-8, 'xtol': 1e-8})  # tol=1e-8) #
#     fit.x = np.abs(fit.x)  # deals with negative parameters for log
#     return fit.x
#
#
# def fitDegradationEquation_RFgreaterThanOne(DDD_vs_remaining_factor, parameters=None):
#     """
#     Given the displacement damage dose curve and minimizing using the Nelder-Mead method to arrive at the best fit for A, C, and D_x for the DDD double degradation equation.  The form of the degradation equation used for the fit replaces the 1 with the fitting parameter A.  This is done for data that has a remaining factor above 1 after radiation.  We do not recommend this method as a remaining factor above indicates a measurement or calibration error.  It does however provide a way to fit otherwise unfittable data.  The data is fit by minimizing the error function to the root mean squared error (RMSE).  This is done becuase it gives lower chi squared values and r squared values closer to 1.  Also, you can be pretty far off on your starting parameters and still arrive at a good fit.  RMSE allows for better fits as there is usually error in the fluence and remaining factor and provides the best fit to expected values
#
#     Args:
#         DDD_vs_remaining_factor (ndarray) : 2d numpy array where column 0 is the displacement damage dose and column 1 is the remaining factor
#         parameters (ndarray): 1d numpy array of starting values for the fit. The array should be 3 elements where the first element is the A parameter in the degradation equation, C is the second parameter, and D_x is the third paramter. If no parameters are entered the fit defaults to using 0.9, 2e-1, 2e10 for A, C, and D_x respectively.  These parameters are found to fit most data without any modification
#
#     Returns:
#             1d numpy array of 3 elements where the first element is the A parameter, the second is C, and the third parameter is D_x.
#
#     """
#     if parameters == None:
#         parameters = [0.9, 2e-1, 2e10]
#
#     if DDD_vs_remaining_factor[0, 0] != 0:
#         D_d = np.vstack(([0, 1], DDD_vs_remaining_factor))  # make sure we fit to 0 and 1
#     else:
#         D_d = DDD_vs_remaining_factor
#
#     fit = sp.optimize.minimize(error_function_degradation_equation_RFgreater_than_one(), parameters,
#                                args=(D_d[:, 0], D_d[:, 1]), method='Nelder-Mead',
#                                options={'maxiter': 10000, 'maxfev': 10000, 'ftol': 1e-6, 'xtol': 1e-6})  # tol=1e-8) #
#     fit.x = np.abs(fit.x)  # deals with negative parameters for log
#     return fit.x
#
#
def fit_nValue_effectiveDDD(particle_energy, fluence, remaining_factor, NIEL, energy_to_normalize, fit_type='single'):
    """
    Give the particle energies, fluences, and remaining factors of an irradiated cell this function fits the data to find the n_value that best collapses the displacement damage dose curves.  The function minimizes by fitting the displacement damage dose curves using the degradation equation and n value to find the best bit using r squared as the fitting parameter.  I have found this to provide better fits than those in the SCREAM model fits

    Args:
        particle_energy (ndarray of float) : 1d numpy array of particle energies or single float value of particle energies
        fluence (ndarray) : 1d numpy array of fluences.  The fluences are converted to displacement damage dose using the NIEL curve
        remaining_factor (ndarray) : 1d numpy array of remaining factors for each fluence
        NIEL (ndarray) : 2d numpy array where column 0 is the particle energy and column 1 is the Non Ionizing Energy Loss of the material
        energy_to_normalize (float): particle energy to normalize the effective displacement damage dose curve. This value is typically 1 for 1MeV electrons

    Returns:
        The n value that best collapses the displacement damage dose curves

    """
    fit = sp.optimize.fminbound(error_function_nValue_for_effectiveDDD, 0, 4,
                                args=(particle_energy, fluence, remaining_factor, NIEL, energy_to_normalize, fit_type), xtol=1e-8)
    return fit


# def fit_degradation_eq_leastsquare(DDD_vs_remaining_factor):
#     """
#     Given the displacement damage dose curve and minimizing using sum of least squares to arrive at the best fit for C and D_x for the DDD degradation equation.
#
#     Args:
#         DDD_vs_remaining_factor (ndarray) : 2d numpy array where column 0 is the displacement damage dose and column 1 is the remaining factor
#
#     Returns:
#         1d numpy array of 2 elements where the first element is the C parameter and the second is D_x
#
#     """
#     # fits a little better than my personal fit
#     if DDD_vs_remaining_factor[0, 0] != 0:
#         D_d = np.vstack(([0, 1], DDD_vs_remaining_factor))  # make sure we fit to 0 and 1
#     else:
#         D_d = DDD_vs_remaining_factor
#
#     popt, pcov = sp.optimize.curve_fit(degradation_equation, D_d[:, 0], D_d[:, 1])
#     return popt
#
# def fit_double_degradation_eq_leastsquare(DDD_vs_remaining_factor):
#     """
#     Given the displacement damage dose curve and minimizing using sum of least squares to arrive at the best fit for C1, D_x1, C2, D_x2 for the DDD double degradation equation.  The data is fit by minimizing the error function to the root mean squared error (RMSE).  This is done becuase it gives lower chi squared values and r squared values closer to 1.  Also, you can be pretty far off on your starting parameters and still arrive at a good fit.  RMSE allows for better fits as there is usually error in the fluence and remaining factor and provides the best fit to expected values
#
#     Args:
#         DDD_vs_remaining_factor (ndarray) : 2d numpy array where column 0 is the displacement damage dose and column 1 is the remaining factor.
#
#     Returns:
#         1d numpy array of 2 elements where the first element is the C1 parameter in the degradation equation, D_x1 is the second parameter, C2 is the third parameter, and D_x2 is the last paramenter.
#     """
#     # fits a little better than my personal fit
#     if DDD_vs_remaining_factor[0, 0] != 0:
#         D_d = np.vstack(([0, 1], DDD_vs_remaining_factor))  # make sure we fit to 0 and 1
#     else:
#         D_d = DDD_vs_remaining_factor
#
#     popt, pcov = sp.optimize.curve_fit(double_degradation_equation, D_d[:, 0], D_d[:, 1])
#     return popt


def getTotalDDD(energy_vs_ddd, cutoff=0):
    """
    Takes a given displacement damage dose particle spectrum and calculates the total displacement damage dose expected for a given environment

    Args:
        energy_vs_ddd (ndarray) : 2d numpy array where column 0 is the particle energy and column 1 is the displacment damage does per particle energy

    Returns:
        The total displacement damage does in a given particle environment

    """
    #just a wrapper for scipy's trapezoidal integration
    energy_vs_ddd = energy_vs_ddd[energy_vs_ddd[:,0]>=cutoff]
    totalDDD = sp.integrate.trapz(energy_vs_ddd[:, 1], energy_vs_ddd[:, 0])
    return totalDDD


def get_energy_vs_ddd(particle_spectrum, NIEL, energy_to_normalize=None, n_value=None):
    """
    Takes a give particle spectrum for an environment and calculates the displacement damage dose for that environment using the non-ionizing energy loss curve.

    Args:
        particle_spectrum (ndarray) : 2d numpy array where column 0 is the particle energy and column 1 is the partcle fluence
        NIEL (ndarray) : 2d numpy array where column 0 is the particle energy and column 1 is the non-ionizing energy loss
        energy_to_normalize (float): The energy to relate all displacement damage dose curves to.  For electrons this is typically 1
        n_value (float) : Empirically determined fitting parameter to collapse displacement damage dose curves that do not collapse. It is 1 for protons

    Returns:
        2d numpy array where column 0 is the particle energy and column 1 is the displacement damage dose.

    """
    if energy_to_normalize == None:
        normalizedNIEL = 1
    else:
        normalizedNIEL = lookupNIELValue(energy_to_normalize, NIEL)

    if n_value == None:
        n_value = 1

    ddd = []
    effectiveElectronWeight = 1.0 / (normalizedNIEL ** (n_value - 1))
    # ddd = get_ddd(particle_spectrum[:,0], particle_spectrum[:,1], NIEL)
    for i, differentialFluence in enumerate(particle_spectrum[:, 1]):
        niel = lookupNIELValue(particle_spectrum[i, 0], NIEL)
        if niel < 1e-9: # TODO: Do something smarter here. this is just a stopgap not sure why i did this....need to investigate.  Code could much simpler, no need for the for loop as lookupNIEL can take arrays
            niel = 0
        dddSingle = effectiveElectronWeight * (niel ** n_value) * differentialFluence
        # if particle_spectrum[i,0] == 1:
        #     print('1 MeV DDD calc')
        #     print('fluence:    {}'.format(differentialFluence))
        #     print('niel:    {}'.format(niel))
        #     print('ddd: {}'.format(dddSingle))

        ddd.append(dddSingle)
    ddd = np.array(ddd)
    energy_vs_ddd = np.zeros((np.shape(particle_spectrum)))
    energy_vs_ddd[:, 0] = particle_spectrum[:, 0]
    energy_vs_ddd[:, 1] = ddd
    return energy_vs_ddd


def get_cumulative_ddd(particle_spectrum, NIEL, energy_to_normalize=None, n_value=None, norm=True):
    """
    Return the cumulative displacement damage dose spectrum for a give environment particle spectrum using the non-ionizing energy loss curve

    Args:
        particle_spectrum (ndarray) : 2d numpy array where column 0 is the particle energy and column 1 is the partcle fluence
        NIEL (ndarray) : 2d numpy array where column 0 is the particle energy and column 1 is the non-ionizing energy loss
        energy_to_normalize (float: The energy to relate all displacement damage dose curves to.  For electrons this is typically 1
        n_value (float) : Empirically determined fitting parameter to collapse displacement damage dose curves that do not collapse. It is 1 for protons
        norm (bool): Optional argument to normalize the cumulative ddisplacement damage dose to 1 if set to True

    Returns:
        2d numpy array where column 0 is the particle energy and column 1 is the cumulative displacement damage dose

    """
    energy_vs_ddd = get_energy_vs_ddd(particle_spectrum, NIEL, energy_to_normalize, n_value)
    cumulative_ddd = np.zeros(np.shape(particle_spectrum))[:-1]
    cumulative_ddd[:, 0] = particle_spectrum[:-1, 0]
    cumulative_ddd[:, 1] = sp.integrate.cumtrapz(energy_vs_ddd[:, 1], energy_vs_ddd[:, 0])
    if norm:
        cumulative_ddd[:, 1] = cumulative_ddd[:, 1] / np.max(cumulative_ddd[:, 1])
    return cumulative_ddd


def interpolateNIELvsTd(newTd, NIELs_vs_Td, Tds, particle_energies, kind='linear'):
    """
    By using 2d interpolation of a table of NIEL values for various Tds, this function will return a NIEL curve at a specific Td using 2d linear interpolation.  The appropriate method would simple be to caluculate the NIEL curve instead of doing linear interpolation. This funtion will be replaced with the calculation in the future

    Args:
        newTd (float): Threshold displacement energy to derive a new non-ionizing energy loss curve
        NIELs_vs_Td (ndarray): 2d numpy array of NIEL values only where the x values or colums are the threshold displacement energies used to calculate the NIEL curves and the y values or rows are the particle energies of the NIEL curves
        Tds (ndarray): 1d numpy array of the threshold displacement energies of the NIEL values in the NIELs_vs_Td data set
        particle_energies (ndarray): 1d numpy array of the particle energies of the NIEL values in the NEIL_vs_Td data set
        kind (str) : The type if interpolation to use. The default is 'linear".  Refer to documentation on scipy interp2d for other interpolation methods

    Returns:
        2d numpy arrary where column 0 is the particle energy and column 1 is the NIEL curve using the user defined threshold displacement energy

    """
    interpolateFunction = sp.interpolate.interp2d(Tds, particle_energies, NIELs_vs_Td, kind=kind)
    newNIEL = interpolateFunction(newTd, particle_energies)
    newNIEL = newNIEL.flatten()
    newTdNIEL = np.zeros((len(newNIEL),2))
    newTdNIEL[:,0] = np.copy(particle_energies)
    newTdNIEL[:,1] = np.copy(newNIEL)
    return newTdNIEL

def errorFunctionBestThresholdDisplacement(T_d, NIELs_vs_Td, Tds, particle_energies_NIEL, particle_energies_qual_data, fluence, remaining_factor):
    NIELcurve = interpolateNIELvsTd(T_d, NIELs_vs_Td, Tds, particle_energies_NIEL)
    ddd_vs_remaining_factor = get_ddd_vs_remaining_factor(particle_energies_qual_data, fluence, remaining_factor, NIELcurve)
    fit = eq.fit_degradation_equation(ddd_vs_remaining_factor)
    newRemainingFactors = eq.degradation_equation(ddd_vs_remaining_factor[:,0], *fit)
    goodnessOfFit = analyze_goodness_of_fit(ddd_vs_remaining_factor[:, 1], newRemainingFactors)
    error = 1-goodnessOfFit.r2 # TODO: think about use chi_squared
    return error


def errorFunctionBestThresholdDisplacement_doubleDegradation(T_d, NIELs_vs_Td, Tds, particle_energies_NIEL, particle_energies_qual_data, fluence, remaining_factor):
    NIELcurve = interpolateNIELvsTd(T_d, NIELs_vs_Td, Tds, particle_energies_NIEL)
    ddd_vs_remaining_factor = get_ddd_vs_remaining_factor(particle_energies_qual_data, fluence, remaining_factor, NIELcurve)
    fit = eq.fit_double_degradation_equation(ddd_vs_remaining_factor)
    newRemainingFactors = eq.double_degradation_equation(ddd_vs_remaining_factor[:,0], *fit)
    goodnessOfFit = analyze_goodness_of_fit(ddd_vs_remaining_factor[:, 1], newRemainingFactors)
    error = 1-goodnessOfFit.r2 # TODO: think about use chi_squared
    return error


def fitToFindThresholdDisplacement(NIELs_vs_Td, Tds, particle_energies_NIEL, particle_energies_qual_data, fluence, remaining_factor):
    """
    Fit to find the best NIEL curve to fit displacement damage dose data that does not collapse using the single degradation equation.  The NIEL curves are selected based on the threshold displacement energy (Td) used to derive the NIEL curve.  This function takes a 2d table where x is Td and y is the particle energy to do interpolation to determine the best NIEL curve on Td.  Once the appropriate NIEL Td curve has been determined the actual NIEL curve using the Td derived from this fit should be used to check if it is correct. This is because 2d linear interpolation is used to derive the appropriate NIEL curve.  So far it looks to fit pretty well

    Args:
        NIELs_vs_Td: 2d numpy array of NIEL values only where the x values or columns are the threshold displacement energies used to calculate the NIEL curves and the y values or rows are the particle energies of the NIEL curves
        Tds: 1d numpy array of the threshold displacement energies of the NIEL values in the NIELs_vs_Td data set
        particle_energies_NIEL: 1d numpy array of the particle energies of the NIEL values in the NEIL_vs_Td data set
        particle_energies_qual_data: 1d numpy array of the particle energies of the radiation data to be collapsed using displacement damage dose method
        fluence: 1d numpy array of of the fluences of the radiation data to be collapsed using displacement damage dose method
        remaining_factor: of the remaining factors of the radiation data to be collapsed using displacement damage dose method

    Returns:
        Threshold displament energy of the NIEL curve that collapse the radiation data using a n value of 1

    """
    minBound = 0
    maxBound = 40
    fit = sp.optimize.fminbound(errorFunctionBestThresholdDisplacement, minBound, maxBound, args= (NIELs_vs_Td, Tds, particle_energies_NIEL, particle_energies_qual_data, fluence, remaining_factor), xtol=1e-10)
    return fit

def fitToFindThresholdDisplacement_double_degradation(NIELs_vs_Td, Tds, particle_energies_NIEL, particle_energies_qual_data, fluence, remaining_factor):
    """
    Fit to find the best NIEL curve to fit displacement damage dose data that does not collapse using the double degradation equation.  The NIEL curves are selected based on the threshold displacement energy (Td) used to derive the NIEL curve.  This function takes a 2d table where x is Td and y is the particle energy to do interpolation to determine the best NIEL curve on Td.  Once the appropriate NIEL Td curve has been determined the actual NIEL curve using the Td derived from this fit should be used to check if it is correct. This is because 2d linear interpolation is used to derive the appropriate NIEL curve.  So far it looks to fit pretty well

    Args:
        NIELs_vs_Td: 2d numpy array of NIEL values only where the x values or colums are the threshold displacement energies used to calculate the NIEL curves and the y values or rows are the particle energies of the NIEL curves
        Tds: 1d numpy array of the threshold displacement energies of the NIEL values in the NIELs_vs_Td data set
        particle_energies_NIEL: 1d numpy array of the particle energies of the NIEL values in the NEIL_vs_Td data set
        particle_energies_qual_data: 1d numpy array of the particle energies of the radiation data to be collapsed using displacement damage dose method
        fluence: 1d numpy array of of the fluences of the radiation data to be collapsed using displacement damage dose method
        remaining_factor: of the remaining factors of the radiation data to be collapsed using displacement damage dose method

    Returns:
        Threshold displament energy of the NIEL curve that collapse the radiation data using a n value of 1

    """
    minBound = 0
    maxBound = 40
    fit = sp.optimize.fminbound(errorFunctionBestThresholdDisplacement_doubleDegradation, minBound, maxBound, args= (NIELs_vs_Td, Tds, particle_energies_NIEL, particle_energies_qual_data, fluence, remaining_factor), xtol=1e-8)
    return fit

def get_ddd_vs_remaining_factor_by_particle_energy(qual_data, NIEL, energy_to_normalize=None, n_value=None):
    """
    Organizes a table of radiation qualification data into groups based on particle energy.  This function is useful for organizing radiation qualification data tables for plotting and analysis

    Args:
        qual_data (ndarray): 2d numpy array where column 0 is hte particle energy for each fluence, column 1 is the fluence, and column 2 is the remaining factor
        NIEL (ndarray) : 2d numpy array where column 0 is the particle energy and column 1 is the non-ionizing energy loss
        energy_to_normalize (float: The energy to relate all displacement damage dose curves to.  For electrons this is typically 1
        n_value (float) : Empirically determined fitting parameter to collapse displacement damage dose curves that do not collapse. It is 1 for protons

    Returns:
        A list of 2d numpy arrays where each element in the list represents one particle energy and and each array in the list is composed of three columns where column 1 is the particle energy, column 2 is the fluence, and column 3 is the remaining factor

    """
    grouped_by_energy = _groupByEnergy(qual_data[:,0], qual_data)
    ddd_vs_remaining_factor_by_particle_energy = []
    for energy in grouped_by_energy:
        ddd_vs_remaining_factor_by_particle_energy.append(get_ddd_vs_remaining_factor(energy[:,0], energy[:,1], energy[:,2], NIEL, energy_to_normalize, n_value))
    return ddd_vs_remaining_factor_by_particle_energy


def _groupByEnergy(grouped_by, data_to_group):
    groupedByEnergy = []
    indicesOfDuplicates = []
    groups = np.unique(grouped_by)
    for group in groups:
        indicesOfDuplicates.append(list(np.where(grouped_by == group)))
    for index in indicesOfDuplicates:
        groupedByEnergy.append(data_to_group[np.min(index):np.max(index) + 1])
    return groupedByEnergy


# proton to electrons to protons back to electrons
def convertDDDelectronsToDDDprotons(D_px, D_ex, C_e, C_p, D_e, A_p=1, A_e=1):
    D_eTOp = D_px * (((1 + (D_e / D_ex)) ** (C_e / C_p)) - 1)
    return D_eTOp


def convertDDDprotonsToDDDelectrons(D_px, D_ex, C_e, C_p, D_p, A_p=1, A_e=1):
    D_pTOe = D_ex * (((1 + (D_p / D_px)) ** (C_p / C_e)) - 1)
    return D_pTOe

def adjusted_niel(energy_0, niel_curve, range_table, thickness_cm):
    intp_range_table = sp.interpolate.interp1d(range_table[:,0], range_table[:,1], fill_value='extrapolate')
    intp_energy_table = sp.interpolate.interp1d(range_table[:,1], range_table[:,0], fill_value='extrapolate')
    adj_niel = []
    valid_niel_energies = []
    # probably can do this with list comprehension
    for energy in energy_0:
        # if (niel_energy > np.min(range_table[:,0])) and (niel_energy < np.max(range_table[:,0])):
        range_niel_energy = intp_range_table(energy)
        if range_niel_energy <= thickness_cm:
            valid_niel_energies.append(energy)
            # if energy >= 0.05:
                # print()
            E_n = niel_curve[0,0]
            # E_n = 0.0002
            j_n = np.log10(E_n / energy) / np.log10(0.999)
            j = np.arange(0,j_n, 1)
            E_j = (0.999 ** j) * energy
            E_j = E_j[E_j >= np.min(range_table[:,0])]
            # E_j = E_j[E_j <= np.max(range_table)]
            # adj_niel.append(np.sum(lookupNIELValue(E_j[:-1], niel_curve) * (np.diff(intp_range_table(E_j[::-1]))[::-1])) / thickness_cm)
            # adj_niel.append(np.sum(lookupNIELValue(E_j, niel_curve) * np.abs(np.diff(intp_range_table(E_j),append=0))) / thickness_cm)
            depth = np.cumsum(np.abs(np.diff(intp_range_table(E_j),append=0)))
            adj_niel.append((sp.integrate.trapz(lookupNIELValue(E_j, niel_curve)*5.32, depth)/thickness_cm)/5.32)
            # adj_niel.append(np.sum(lookupNIELValue(E_j[:-1], niel_curve) * np.abs(np.diff(intp_range_table(E_j))) / thickness_cm))

        elif (range_niel_energy > thickness_cm):
            valid_niel_energies.append(energy)
            E_n = intp_energy_table(intp_range_table(energy) - thickness_cm)
            j_n = np.log10(E_n / energy) / np.log10(0.99)
            # j_n = np.log10(E_n) / np.log10(0.99)
            j = np.linspace(0,j_n)
            E_j = (0.99 ** j) * energy
            E_j = E_j[E_j >= np.min(range_table[:,0])]
            # E_j = E_j[E_j <= np.max(range_table)]
            adj_niel.append(np.sum(lookupNIELValue(E_j[:-1], niel_curve) * np.abs(np.diff(intp_range_table(E_j)))) / thickness_cm)

    adj_niel_2d = np.vstack((valid_niel_energies, adj_niel)).T
    # print(adj_niel_2d)
    return adj_niel_2d

        # j_n = np.array(j_n)
        # E_j = (0.99 ** j_n) * niel_energy
        # range_E_j = intp_range_table(E_j)
        # adj_niel = np.sum(lookupNIELValue(E_j)*np.diff(range_E_j))/thickness_cm
        # j_n = []
        # print(niel_energy, E_n, j_n, 0.99**j_n * niel_energy)

        # print(E_j)
    # return E_j

def convert_ddd_to_fluence(ddd, particle_energy, niel_curve):
    niel = lookupNIELValue(particle_energy, niel_curve)
    fluence = ddd/niel
    return fluence
