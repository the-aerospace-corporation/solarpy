import numpy as np

class analyze_goodness_of_fit:
    def __init__(self, observedValues, predictedValues):
        self.observedValues = observedValues
        self.predictedValues = predictedValues
        self.r2 = self._rSquared()
        self.residuals = self._residuals()
        self.sumOfSquaresofResiduals = self._sumOfSquaresofResiduals()

    def _rSquared(self):
        r2 = 1 - np.abs(self._sumOfSquaresofResiduals() / self._totalSumofSquares())
        return r2

    def _totalSumofSquares(self):
        SS_total = (np.abs(self.observedValues - self._meanOfObservedValues()) ** 2).sum()
        return SS_total

    def _regressionSumOfSquares(self):
        regressionSumOfSquares = (np.abs(self.observedValues - self._meanOfObservedValues()) ** 2).sum()
        return regressionSumOfSquares

    def _sumOfSquaresofResiduals(self):
        sumOfSquaresOfResiduals = (np.abs(self.observedValues - self.predictedValues) ** 2).sum()
        return sumOfSquaresOfResiduals

    def _residuals(self):
        residuals = self.observedValues - self.predictedValues
        return residuals

    def _meanOfObservedValues(self):
        meanOfObservedValues = np.mean(self.observedValues)
        return meanOfObservedValues
