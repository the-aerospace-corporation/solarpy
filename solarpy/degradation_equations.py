import numpy as np
import scipy as sp
from scipy import interpolate, optimize
import matplotlib.pyplot as plt

from solarpy.analyze_goodness_of_fit import analyze_goodness_of_fit


class singleDegradationEquation(object):
    def __init__(self, fluence_vs_remaining_factor):
        self.fluence_vs_remaining_factor = fluence_vs_remaining_factor
        self.coefficients = None

    def fit(self, parameters=None):
        if parameters is None:
            parameters = [2e-2, 2e8]
        coefficients = fit_degradation_equation(self.fluence_vs_remaining_factor, parameters=parameters)
        self.r2 = analyze_goodness_of_fit(self.fluence_vs_remaining_factor[:,1], degradation_equation(self.fluence_vs_remaining_factor[:,0],*coefficients)).r2
        self.coefficients = coefficients
        return coefficients

    def getFluence(self, remaining_factor):
        if self.coefficients is None:
            self.coefficients = self.fit()
        fluence = get_fluence(remaining_factor, *self.coefficients)
        return fluence

    def getRemainingFactor(self, fluence):
        if self.coefficients is None:
            self.coefficients = self.fit()
        remaining_factor = degradation_equation(fluence, *self.coefficients)
        return remaining_factor



class doubleDegradationEquation(object):
    def __init__(self, fluence_vs_remaining_factor, zero_fit='n'):
        self.fluence_vs_remaining_factor = fluence_vs_remaining_factor
        self.coefficients = None
        if zero_fit == 'y':
            self.fluence_vs_remaining_factor = np.vstack(([0,1], self.fluence_vs_remaining_factor))

    def fit(self, parameters=None):
        if parameters is None:
            parameters = [5e-3, 5e10, 5e-3, 5e-3]
        coefficients = fit_double_degradation_equation(self.fluence_vs_remaining_factor, parameters=parameters)
        self.r2 = analyze_goodness_of_fit(self.fluence_vs_remaining_factor[:,1], double_degradation_equation(self.fluence_vs_remaining_factor[:,0],*coefficients)).r2
        self.coefficients = coefficients
        return coefficients

    def getFluence(self, remaining_factor):
        if self.coefficients is None:
            self.coefficients = self.fit()
        fluence = get_fluence_ddd_double_degradation(remaining_factor, *self.coefficients)
        return fluence

    def getRemainingFactor(self, fluence):
        if self.coefficients is None:
            self.coefficients = self.fit()
        remaining_factor = double_degradation_equation(fluence, *self.coefficients)
        return remaining_factor


def degradation_equation(fluence_ddd, C, phi_D_x):
    """
    Degradation equation as defined by the Solar Cell and GaAs Solar Cell Radiation Handbooks (1 - C * np.log10(1 + (fluence / phi_x))). The same equation is used for displacement damage dose curves

    Args:
        fluence_ddd (ndarray): 1d numpy array of fluences or displacement damage dose
        C (float): Constant that is determined through fitting
        phi_D_x(float): Fluence or displacement damage dose constant that is determined through fitting

    Returns:
        1d array of calculated remaining factors
    """
    phi_D_x = np.abs(phi_D_x)  # prevents negative numbers into log
    remainingFactor = 1 - C * np.log10(1 + (fluence_ddd / phi_D_x))
    return remainingFactor


def double_degradation_equation(fluence_ddd, C1, phi_D_x1, C2, phi_D_x2):
    """
    Combines two of the single degradation equation (1 - (C1 * np.log10(1 + (fluence / phi_x1))) - (C2 * np.log10(1 + (fluence / phi_x2)))) in an effort to alleviate the assumption that there are two components degrading.  Can fit more complex solar cell degradation better, but may not always be the aright approach as as it is an assumption that there are two mechanisms behaving similarly

    Args:
        fluence_ddd (ndarray): 1d numpy array of fluences or displacement damage dose
        C1 (float) : Constant that is determined through fitting for first degradation equation
        phi_D_x1 (float) : Fluence or displacement damage dose constant that is determined through fitting
        C2 (float) : Constant that is determined through fitting for first degradation equation
        phi_D_x2 (float) : Fluence or displacement damage dose constant that is determined through fitting

    Returns:
        1d array of calculated remaining factors

    """
    phi_D_x1 = np.abs(phi_D_x1)  # prevents negative numbers into log
    phi_D_x2 = np.abs(phi_D_x2)
    remainingFactor = 1 - (C1 * np.log10(1 + (fluence_ddd / phi_D_x1))) - (C2 * np.log10(1 + (fluence_ddd / phi_D_x2)))
    return remainingFactor


def degradation_equation_rf_greater_than_one(fluence_ddd, A, C, phi_D_x):
    """
    Degradation equation for determining remaining factor as a function of displacement damage dose.  The 1 in the degradation equation is replaced with a fitting parameter (A) to account for values above one.  This should not be used because it assumes that the remaining factors can be greater than beginning of life. It is in included as a reference to previous literature.

    Args:
        fluence_ddd (ndarray): 1d numpy array of fluences or displacement damage dose
        A (float) : Constant that is determined from fitting to account for remaining factors greater than 1
        C (float) : Constant that is determined through fitting
        phi_D_x(float): Fluence or displacement damage dose constant that is determined through fitting

    Returns:
        1d array of calculated remaining factors

    """
    phi_D_x = np.abs(phi_D_x)  # prevents negative numbers into log
    remainingFactor = A - C * np.log10(1 + (fluence_ddd / phi_D_x))
    return remainingFactor

def double_degradation_equation_rf_greater_than_one(fluence_ddd, A, C1, phi_D_x1, C2, phi_D_x2):
    """
    Combines two of the single degradation equation (A - (C1 * np.log10(1 + (fluence / phi_x1))) - (C2 * np.log10(1 + (fluence / phi_x2)))) in an effort to alleviate the assumption that there are two components degrading.  Can fit more complex solar cell degradation better, but may not always be the aright approach as as it is an assumption that there are two mechanisms behaving similarly

    Args:
        fluence_ddd (ndarray): 1d numpy array of fluences or displacement damage dose
        C1 (float) : Constant that is determined through fitting for first degradation equation
        phi_D_x1 (float) : Fluence or displacement damage dose constant that is determined through fitting
        C2 (float) : Constant that is determined through fitting for first degradation equation
        phi_D_x2 (float) : Fluence or displacement damage dose constant that is determined through fitting

    Returns:
        1d array of calculated remaining factors

    """
    phi_D_x1 = np.abs(phi_D_x1)  # prevents negative numbers into log
    phi_D_x2 = np.abs(phi_D_x2)
    remainingFactor = A - (C1 * np.log10(1 + (fluence_ddd / phi_D_x1))) - (C2 * np.log10(1 + (fluence_ddd / phi_D_x2)))
    return remainingFactor

def errorFunctionDegradationEquation(p, fluence_ddd, remaining_factor):
    number_of_y_points = len(fluence_ddd)
    RF_new = (((np.abs(remaining_factor - degradation_equation(fluence_ddd, *p)) / remaining_factor) * 100) ** 2)
    error = np.sqrt(RF_new.sum() / number_of_y_points)
    return error


def errorFunctionDoubleDegradationEquation(p, fluence, remainingFactor):
    number_of_y_points = len(fluence)
    RF_new = (((np.abs(remainingFactor - double_degradation_equation(fluence, *p)) / remainingFactor) * 100) ** 2)
    error = np.sqrt(RF_new.sum() / number_of_y_points)
    return error


def error_function_degradation_equation_rf_greater_than_one(p, fluence_ddd, remaining_factor):
    number_of_y_points = len(fluence_ddd)
    RF_new = (((np.abs(
        remaining_factor - degradation_equation_rf_greater_than_one(fluence_ddd, *p)) / remaining_factor) * 100) ** 2)
    error = np.sqrt(RF_new.sum() / number_of_y_points)
    return error


def fit_degradation_equation(fluence_or_ddd_vs_remaining_factor, parameters=None):
    """
    Given the displacement damage dose or fluence vs remaining factor, the data is fit using the Nelder-Mead method to arrive at the best fit for C and phi_D_x for the degradation equation.  The data is fit by minimizing the error function to the root mean squared error (RMSE).  This is done because it gives lower chi squared values and r squared values closer to 1.  Also, you can be pretty far off on your starting parameters and still arrive at a good fit.  RMSE allows for better fits as there is usually error in the fluence and remaining factors.

    Args:
        fluence_or_ddd_vs_remaining_factor (ndarray) : 2d numpy array where column 0 is the displacement damage dose or fluence and column 1 is the remaining factor:
        parameters (ndarray): 1d numpy array of starting values for the fit. The array should be 2 elements where the first element is the C parameter in the degradation equation and phi_D_x is the second parameter. If no parameters are entered the fit defaults to using 2e-2 for C and 2e8 for phi_D_x.  These parameters are found to fit most data without any modification

    Returns:
        1d numpy array of 2 elements where the first element is the C parameter and the second is phi_D_x

    """
    if parameters is None:
        parameters = [2e-2, 2e8]

    if fluence_or_ddd_vs_remaining_factor[0, 0] != 0:
        fluence_or_ddd_vs_remaining_factor = np.vstack(
            ([0, 1], fluence_or_ddd_vs_remaining_factor))  # make sure we fit to 0 and 1

    fit = sp.optimize.minimize(errorFunctionDegradationEquation, parameters, args=(
        fluence_or_ddd_vs_remaining_factor[:, 0], fluence_or_ddd_vs_remaining_factor[:, 1]), method='Nelder-Mead',
                               tol=1e-9)  # options = {'maxiter':10000, 'maxfev':10000, 'ftol': 1e-6, 'xtol': 1e-6})
    fit = np.abs(fit.x)  # deals with negative parameters for log
    return fit


def fit_double_degradation_equation(fluence_or_ddd_vs_remaining_factor, parameters=None):
    """
    Given the displacement damage dose or fluence vs remaining factor, the data is fit using the Nelder-Mead method to arrive at the best fit for C1, phi_D_x1, C2, phi_D_x2 for the double degradation equation.  The data is fit by minimizing the error function to the root mean squared error (RMSE).  This is done becuase it gives lower chi squared values and r squared values closer to 1.  Also, you can be pretty far off on your starting parameters and still arrive at a good fit.  RMSE allows for better fits as there is usually error in the fluence and remaining factor and provides the best fit to expected values.

    Args:
        fluence_or_ddd_vs_remaining_factor (ndarray) : 2d numpy array where column 0 is the displacement damage dose or fluence and column 1 is the remaining factor:
        parameters (ndarray): 1d numpy array of starting values for the fit. The array should be 4 elements where the first element is the C1 parameter in the degradation equation, D_x1 is the second parameter, C2 is the third parameter, and D_x2 is the last parameter.  If no parameters are entered the fit defaults to using 6e-2, 6e8, 2e-4, 4e-10 for C1, D_x1, C2, and D_x2 respectively.  These parameters are found to fit most data without any modification

    Returns:
        1d numpy array of 4 elements where the first element is the C1 parameter in the degradation equation, phi_D_x1 is the second parameter, C2 is the third parameter, and phi_D_x2 is the last parameter.

    """
    # Best guesses for JPL EQFLUX Fits[1.8e-6, 1e-2, 1.7e-1, 5e9][5.0e-6, 4e8, 2.0e-1, 1e10][1.0e-3, 4e9, 1e-1, 8e9][2.1e-4, 1e9, 2.0e-1, 4e9]
    if parameters is None:
        parameters = [6e-2, 6e9, 2e-5, 4e-10]

    if fluence_or_ddd_vs_remaining_factor[0, 0] != 0:
        fluence_or_ddd_vs_remaining_factor = np.vstack(
            ([0, 1], fluence_or_ddd_vs_remaining_factor))  # make sure we fit to 0 and 1

    fit = sp.optimize.minimize(errorFunctionDoubleDegradationEquation, parameters, args=(
        fluence_or_ddd_vs_remaining_factor[:, 0], fluence_or_ddd_vs_remaining_factor[:, 1]), method='Nelder-Mead',
                               tol=1e-10)  # , options = {'maxiter':10000, 'maxfev':10000, 'ftol': 1e-16, 'xtol': 1e-16})
    fit = np.abs(fit.x)  # deals with negative parameters for log
    return fit


def fit_degradation_equation_rf_greater_than_1(fluence_or_ddd_vs_remaining_factor, parameters=None):
    """
    Given the displacement damage dose or fluence vs remaining factor, the data is fit using the Nelder-Mead method to arrive at the best fit for A, C, and phi_D_x for the DDD double degradation equation.  The form of the degradation equation used for the fit replaces the 1 with the fitting parameter A.  This is done for data that has a remaining factor above 1 after radiation.  We do not recommend this method as a remaining factor above indicates a measurement or calibration error.  It does however provide a way to fit otherwise unfittable data.  The data is fit by minimizing the error function to the root mean squared error (RMSE).  This is done becuase it gives lower chi squared values and r squared values closer to 1.  Also, you can be pretty far off on your starting parameters and still arrive at a good fit.  RMSE allows for better fits as there is usually error in the fluence and remaining factor and provides the best fit to expected values

    Args:
        fluence_or_ddd_vs_remaining_factor (ndarray) : 2d numpy array where column 0 is the displacement damage dose or fluence and column 1 is the remaining factor:
        parameters (ndarray): 1d numpy array of starting values for the fit. The array should be 3 elements where the first element is the A parameter in the degradation equation, C is the second parameter, and D_x is the third paramter. If no parameters are entered the fit defaults to using 0.9, 2e-1, 2e10 for A, C, and phi_D_x respectively.  These parameters are found to fit most data without any modification

    Returns:
            1d numpy array of 3 elements where the first element is the A parameter, the second is C, and the third parameter is D_x.

    """
    if parameters == None:
        parameters = [0.9, 2e-1, 2e10]

    if fluence_or_ddd_vs_remaining_factor[0, 0] != 0:
        fluence_or_ddd_vs_remaining_factor = np.vstack(
            ([0, 1], fluence_or_ddd_vs_remaining_factor))  # make sure we fit to 0 and 1

    fit = sp.optimize.minimize(error_function_degradation_equation_rf_greater_than_one, parameters, args=(
        fluence_or_ddd_vs_remaining_factor[:, 0], fluence_or_ddd_vs_remaining_factor[:, 1]), method='Nelder-Mead',
                               options={'maxiter': 10000, 'maxfev': 10000, 'ftol': 1e-6, 'xtol': 1e-6})  # tol=1e-8) #
    fit.x = np.abs(fit.x)  # deals with negative parameters for log
    return fit.x


def fit_degradation_equation_leastsquare(fluence_or_ddd_vs_remaining_factor):
    """
    Given the displacement damage dose or fluence vs remaining factor, the data is fit using the sum of least squares to arrive at the best fit for C and D_x for the DDD degradation equation.

    Args:
        fluence_or_ddd_vs_remaining_factor (ndarray) : 2d numpy array where column 0 is the displacement damage dose or fluence and column 1 is the remaining factor:

    Returns:
        1d numpy array of 2 elements where the first element is the C parameter and the second is phi_D_x

    """
    # fits a little better than my personal fit

    if fluence_or_ddd_vs_remaining_factor[0, 0] != 0:
        fluence_or_ddd_vs_remaining_factor = np.vstack(
            ([0, 1], fluence_or_ddd_vs_remaining_factor))  # make sure we fit to 0 and 1

    popt, pcov = sp.optimize.curve_fit(degradation_equation, fluence_or_ddd_vs_remaining_factor[:, 0],
                                       fluence_or_ddd_vs_remaining_factor[:, 1])
    popt = np.abs(popt)
    return popt


def fit_double_degradation_eq_leastsquare(fluence_or_ddd_vs_remaining_factor):
    """
    Given the displacement damage dose curve and minimizing using sum of least squares to arrive at the best fit for C1, D_x1, C2, D_x2 for the DDD double degradation equation.  The data is fit by minimizing the error function to the root mean squared error (RMSE).  This is done becuase it gives lower chi squared values and r squared values closer to 1.  Also, you can be pretty far off on your starting parameters and still arrive at a good fit.  RMSE allows for better fits as there is usually error in the fluence and remaining factor and provides the best fit to expected values

    Args:
        fluence_or_ddd_vs_remaining_factor (ndarray) : 2d numpy array where column 0 is the displacement damage dose and column 1 is the remaining factor.

    Returns:
        1d numpy array of 2 elements where the first element is the C1 parameter in the degradation equation, D_x1 is the second parameter, C2 is the third parameter, and D_x2 is the last paramenter.
    """
    # fits a little better than my personal fit
    if fluence_or_ddd_vs_remaining_factor[0, 0] != 0:
        fluence_or_ddd_vs_remaining_factor = np.vstack(
            ([0, 1], fluence_or_ddd_vs_remaining_factor))  # make sure we fit to 0 and 1

    popt, pcov = sp.optimize.curve_fit(double_degradation_equation, fluence_or_ddd_vs_remaining_factor[:, 0],
                                       fluence_or_ddd_vs_remaining_factor[:, 1])
    popt = np.abs(popt)
    return popt


def polynomialFit(fluence_or_ddd_vs_remaining_factor, degree=5):
    """
    Given the displacement damage dose or fluence vs remaining factor, the data is fit using a polynomial to arrive at the best fit for C and D_x for the DDD degradation equation.  Wrapper function of numpy's polyfit.

    Args:
        fluence_or_ddd_vs_remaining_factor (ndarray) : 2d numpy array where column 0 is the displacement damage dose or fluence and column 1 is the remaining factor:
        degree (int) : Degree of fitting polynomial

    Returns:
        Polynomial coefficients, highest power first

    """
    if fluence_or_ddd_vs_remaining_factor[0, 0] != 0:
        fluence_or_ddd_vs_remaining_factor = np.vstack(
            ([0, 1], fluence_or_ddd_vs_remaining_factor))  # make sure we fit to 0 and 1

    fit = np.polyfit(np.log(fluence_or_ddd_vs_remaining_factor[:, 0]), fluence_or_ddd_vs_remaining_factor[:, 1], degree)
    return fit


def get_fluence_ddd_poly(remaining_factor, fit, x0=None):
    """
    Calculates the fluence or displacement damage dose using a polynomial fit for a given remaining factor

    Args:
        remaining_factor (ndarray or float) : 1d numpy array or float of the remaining factor
        fit (ndarray) : The polynomial coefficients
        x0 (float) : Best guess at the fluence or displacement damage does that the remaining factor should be near

    Returns:
        1d numpy array or float of the fluence or displacement damage dose that the remaining factor is located at

    """

    if x0 is None:
        x0 = 1e15

    poly = np.poly1d(fit)

    def f(fluence, remainingFactor):
        fluence = np.log(fluence)
        r = remainingFactor - poly(fluence)
        return r

    answer = sp.optimize.fsolve(f, x0=x0, args=(remaining_factor), xtol=1e-8, maxfev=100000000)
    return answer


def get_fluence(remaining_factor, C, phi_D_x):
    """
    Calculates the fluence or displacement damage dose using the degradation equation

    Args:
        remaining_factor (ndarray or float) : 1d numpy array or float of the remaining factor
        C (float): Constant that is determined through fitting
        phi_D_x(float): Fluence or displacement damage dose constant that is determined through fitting

    Returns:
        1d numpy array or float of the fluence or displacement damage dose that the remaining factor is located at

    """
    fluence = phi_D_x * ((10 ** ((1 - remaining_factor) / C)) - 1)
    return fluence


def get_fluence_ddd_double_degradation(remaining_factor, C1, phi_D_x1, C2, phi_D_x2, x0=None):
    """
    Calculates the fluence or displacement damage dose using the double degradation equation

    Args:
        remaining_factor (ndarray or float) : 1d numpy array or float of the remaining factor
        C1 (float) : Constant that is determined through fitting for first degradation equation
        phi_D_x1 (float) : Fluence or displacement damage dose constant that is determined through fitting
        C2 (float) : Constant that is determined through fitting for first degradation equation
        phi_D_x2 (float) : Fluence or displacement damage dose constant that is determined through fitting
        x0 (float) : Best guess at the fluence or displacement damage does that the remaining factor should be near

    Returns:
        1d numpy array or float of the fluence or displacement damage dose that the remaining factor is located at

    """

    if x0 is None:
        fluences = np.logspace(8, 17)
        rf = double_degradation_equation(fluences, C1, phi_D_x1, C2, phi_D_x2)
        fintp = sp.interpolate.interp1d(rf, fluences)
        x0 = fintp(remaining_factor)
        # x0 = 1e11


    def f(fluence, remainingFactor, C1, phi_x1, C2, phi_x2):
        r = remainingFactor - double_degradation_equation(fluence, C1, phi_x1, C2, phi_x2)
        return r

    answer = sp.optimize.fsolve(f, x0=x0, args=(remaining_factor, C1, phi_D_x1, C2, phi_D_x2))[0]
    return answer


def get_remaining_factor_poly(fluence_ddd, fit):
    """
    Calculates the remaining factor from a polynomial fit of the displacement damage dose or fluence.

    Args:
        fluence_ddd (ndarray): 1d numpy array of fluences or displacement damage dose
        fit (ndarray) : The polynomial coefficients

    Returns:
        1d numpy array or float of the remaining factor at a given fluence or displacement damage dose

    """
    poly = np.poly1d(fit)
    remainingFactor = poly(np.log(fluence_ddd))
    return remainingFactor

def plot_fit_checks(qual_data, fit_type=None, fit_parameters=None, normalize_rdc=None, particle_type=None, colors=None, marker=None):
    if fit_type == None:
        fit_type = 'single'

    else:
        fit_type = fit_type.lower()

    if any(isinstance(param, list) for param in fit_parameters):
        multi_param = True
    else:
        multi_param = False

    particle_energy = np.unique(qual_data[:, 0])
    groups = group_by(qual_data, group_by_column_index=0)
    if colors is None:
        colors = ['blue', 'green', 'orange', 'purple', 'pink', 'red', 'cyan', 'midnightblue']
    if marker is None:
        marker = 'o'
    max_x = 0
    min_x = 1e16
    coeff = []
    for i, group in enumerate(groups):
        # plt.plot(group[:, 1], group[:, 2], 'o', color=colors[i], label=str(particle_energy[i]) + 'MeV ' + fit_type)
        plt.plot(group[:, 1], group[:, 2], marker, color=colors[i], label=f'{particle_energy[i]} MeV {particle_type}')

        if multi_param:
            fit_param = fit_parameters[i]
        else:
            fit_param = fit_parameters

        if fit_type == 'single':

            rdc_model = singleDegradationEquation(group[:, [1, 2]])
            rdc_model.fit(parameters=fit_param)
            x = np.logspace(np.log10(np.min(group[:,1])) - 0.2, np.log10(np.max(group[:,1])), 100)
            if i == 0:
                label = 'single-fit'
            else:
                label = None

            plt.plot(x, rdc_model.getRemainingFactor(x), '--', color=colors[i])

        elif fit_type == 'double':
            rdc_model = doubleDegradationEquation(group[:, [1, 2]])
            # rdc_model.fit(parameters=[2e-1, 6e8, 2e-9, 4e-5])
            rdc_model.fit(parameters=fit_param)
            x = np.logspace(np.log10(np.min(group[:,1])) - 0.2, np.log10(np.max(group[:,1])), 100)
            if i == 0:
                label = 'double-fit'
            else:
                label = None

            plt.plot(x, rdc_model.getRemainingFactor(x), '--', color=colors[i])

        if max_x < np.max(x):
            max_x = np.max(x)

        if min_x > np.min(x):
            min_x = np.min(x)

        print(' {} MeV r2: {}, coef: {}'.format(group[0, 0], rdc_model.r2, rdc_model.coefficients))
        coeff.append([group[0, 0]]+list(rdc_model.coefficients))
        # x = np.logspace(np.log10(group[0, 1])-0.2, np.log10(group[-1, 1]), 100)
    # plt.axhline(0.77, linestyle=':', color='black')
    plt.xscale('symlog', linthresh=10**(np.log10(np.min(min_x))-0.1))
    plt.xlim([0, max_x + (0.15*max_x)])
    plt.legend()
    return coeff

def group_by(data, group_by_column_index=0):
    """
    Groups data by a column index
    Args:
        data ():
        group_by_column_index ():

    Returns:

    """
    unique_values = np.unique(data[:,group_by_column_index])
    indicesOfDuplicates = []
    for value in unique_values:
        indicesOfDuplicates.append(list(np.where(data[:,group_by_column_index] == value)[0]))
    grouped_by_energy = []
    for index in indicesOfDuplicates:
        grouped_by_energy.append(data[min(index):max(index) + 1])
    return grouped_by_energy
