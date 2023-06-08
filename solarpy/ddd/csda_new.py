""" degraded_spectra.py
Routines for computing degraded proton or electron spectra for slab geometries
Continuous Slowing Down Approximation (CSDA) only
Paul O'Brien

Note: while CSDA is a very good approximation for protons, it is not
nearly as good for electrons. Electrons scatter significantly.

SCREAM proton and electron range tables provided by Dr. Scott Messenger. SCREAM
proton range data come from SRIM and electron range data from NIST/ESTAR.

NIST/PSTAR proton range tables are also provided.


J(Eout) = integral(0 to pi/2) da J(Ein) * dR/dE(Eout) / (dR/dE(Ein)) sin(a)/2
J(Ein) - incident isotropic differential flux, #/cm^2/MeV
J(Eout) - exiting isotropic differential flux, #/cm^2/MeV
dR/dE(E) - derivative of range w.r.t energy


Initialize a ScreamTable which holds range tables from
ScreamRangeData_protons.csv (default) or ScreamRangeData_electrons.csv
NISTRangeData_protons.csv is also an option (with higher start and end energy limits)
Then initialize a SlabLayer with a material and thickness
Then use the SlabLayer's degraded_spectrum method to compute the degraded
spectrum from incident differential spectrum.

Numerics:

Spectrum is interpolated via power law, log(fluence) vs log(E), except when
fluence has zeros, in which case a semi-log interpolation is used: fluence vs log(E)

Range is interpolated via power law, log(R) vs log(E), and vice versa

Range is differentiated numerically via power law: dR/dE = dlogR/dlogE*R/E
and then interpolated via power law, log(dR/dE) vs log(E)

For energies beyond range table (1 GeV for SCREAM, 10 GeV for NIST),
particle is assumed to traverse any depth of shielding with negligible
energy loss.

Integration is done adaptively using scipy.integrate.quad to relative
default precision of 1e-4.


"""

import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad

_HIGH_ENERGY_WARNED = False  # has high energy pass-thru warning been issued?


class RangeTable(object):
    """
    RangeTable base class
    methods:
        EtoR - energy to range
        RtoE - range to energy
        dRdE - derivative of range w.r.t energy
        density - return density of material
        Rlimits - return range limits as [min,max] in g/cm^2
        Elimits - valid energy limits as [min,max], MeV
    data:
        materials - list of materials

    """

    def __init__(self):
        self.materials = []
        pass  # base class. this must be overloaded in derived class

    def material_check(self, material):
        """
        material_check(material) - check material is known, raise Exception otherwise
        """
        if material not in self.materials:
            raise Exception('Unrecognized material "%s"' % material)
        return True

    def density(self, material):
        """
        density(material) - return material density, g/cm^3
        """
        pass  # base class. this must be overloaded in derived class

    def Elimits(self, material):
        """
        Elimits(material)
        return [min,max] in MeV allowed for material
        """
        pass  # base class. this must be overloaded in derived class

    def Elimit_check(self, material, E):
        """
        Elimit_check(material,E) check E MeV against Elimits for material. Raise exception if out of bounds
        """
        Elimits = self.Elimits(material)
        if np.any(E < Elimits[0]):
            raise Exception('Requested energy %g MeV out of bounds %g-%g' % (np.min(E), Elimits[0], Elimits[1]))
        if np.any(E > Elimits[1]):
            raise Exception('Requested energy %g MeV out of bounds %g-%g' % (np.max(E), Elimits[0], Elimits[1]))
        return True

    def Rlimits(self, material):
        """
        Rlimits(material)
        return [min,max] in g/cm^2 allowed for material
        """
        pass  # base class. this must be overloaded in derived class

    def Rlimit_check(self, material, R):
        """
        Rlimit_check(material,range) check range g/cm^2 against Rlimits for material. Raise exception if out of bounds
        """
        self.material_check(material)
        Rlimits = self.Rlimits(material)
        if np.any(R < Rlimits[0]):
            raise Exception('Requested range %g g/cm^2 out of bounds %g-%g' % (np.min(R), Rlimits[0], Rlimits[1]))
        if np.any(R > Rlimits[1]):
            raise Exception('Requested range %g g/cm^2 out of bounds %g-%g' % (np.max(R), Rlimits[0], Rlimits[1]))
        return True

    def EtoR(self, material, E):
        """
        range = EtoR(material,E) - return range in g/cm^2 given material and energy in MeV
        """
        pass  # base class. this must be overloaded in derived class

    def RtoE(self, material, range):
        """
        E = RtoE(material,range) - return energy in MeV, given material and range in g/cm^2
        """
        pass  # base class. this must be overloaded in derived class

    def dRdE(self, material, E):
        """
        grad = dRdE(material,E) - returns dRdE in g/cm^2/MeV given material and energy
        (derivative of range w.r.t energy)
        """
        pass  # base class. this must be overloaded in derived class


class ScreamTable(RangeTable):
    """ ranges = ScreamTabale(filename=None)
    filename = is the relative or absolute path to the range table file
    default filename is ScreamRangeData_protons.csv (proton range table)
    for electrons, supply ScreamRangeData_electrons.csv
    NIST proton range can be used via NISTRangeData_protons.csv
    data:
        materials - list of materials
    """

    def __init__(self, filename=None):

        super().__init__()

        if filename is None:
            filename = 'ScreamRangeData_protons.csv'
        if not os.path.exists(filename):
            raise Exception('%s does not exist' % filename)
        self.filename = filename

        ranges = {}
        """ranges is a dict
            ranges['Energy (MeV)'] is a list of energies in MeV
            all other keys are materials with values that are themselves dicts
            value dicts have two fields: density and range
            density is a scalar in g/cm^3
            range is alist of g/cm^2 corresponding to entries in Energy list"""

        """density,2.32,2.32,...
           Energy,,,...
           (MeV),Si,Si2,...
           1.000E-04,6.264E-07,"""

        ranges = None
        with open(filename, 'r') as f:
            headers = []
            for i in range(3):
                headers.append(f.readline().strip().split(','))
            densities = np.array([float(x) for x in headers[0][1:]])
            materials = [(a + ' ' + b).strip() for (a, b) in zip(headers[1], headers[2])]
            data = []
            for line in f.readlines():
                line = line.strip();
                if line == ',' * len(line):  # all commas
                    continue
                line = [float(x) for x in line.split(',')]
                data.append(line)
            data = np.array(data)
        ranges = {}
        for i in range(len(materials)):
            if i == 0:
                ranges[materials[0]] = data[:, 0]  # energy
            else:
                rec = {'density': densities[i - 1]}
                rec['range'] = data[:, i]
                ranges[materials[i]] = rec
        ranges['materials'] = materials[1:]  # remove first one, which is energy
        self.ranges = ranges
        self.materials = ranges['materials']
        self._dRdE = {}  # cached gradients, keys are materials

    def Elimits(self, material):
        """
        Elimits(material)
        return [min,max] in MeV allowed for material
        """
        self.material_check(material)
        Elimits = [self.ranges['Energy (MeV)'][0], self.ranges['Energy (MeV)'][-1]]
        return Elimits

    def Rlimits(self, material):
        """
        Rlimits(material)
        return [min,max] in g/cm^2 allowed for material
        """
        self.material_check(material)
        Rlimits = [self.ranges[material]['range'][0], self.ranges[material]['range'][-1]]
        return Rlimits

    def density(self, material):
        """
        density(material) - return material density, g/cm^3
        """
        self.material_check(material)
        return self.ranges[material]['density']

    def EtoR(self, material, E):
        """
        range = EtoR(material,E) - return range in g/cm^2 given material and energy in MeV
        """
        self.material_check(material)
        self.Elimit_check(material, E)
        x = self.ranges['Energy (MeV)']
        y = self.ranges[material]['range']
        f = interp1d(np.log(x), np.log(y), 'linear', copy=False, assume_sorted=True)
        return np.exp(f(np.log(E)))

    def RtoE(self, material, R):
        """
        E = RtoE(material,range) - return energy in MeV, given material and range in g/cm^2
        """
        self.material_check(material)
        self.Rlimit_check(material, R)
        x = self.ranges[material]['range']
        y = self.ranges['Energy (MeV)']
        f = interp1d(np.log(x), np.log(y), 'linear', copy=False, assume_sorted=True)
        return np.exp(f(np.log(R)))

    def dRdE(self, material, E):
        """
        grad = dRdE(material,E) - returns dRdE in g/cm^2/MeV given material and energy
        (derivative of range w.r.t energy)
        """
        self.material_check(material)
        _dRdE = None
        x = self.ranges['Energy (MeV)']
        if material in self._dRdE:  # retrieve precomputed numerical gradient
            _dRdE = self._dRdE[material]
        else:  # compute numerical gradient
            y = self.ranges[material]['range']
            _dRdE = np.gradient(np.log(y), np.log(x)) * y / x  # derivative with power-law interpolation
            assert np.all(_dRdE > 0), 'log-log derivative did not work for material %s' % material
        # linearly interpolate numerical gradient
        f = interp1d(np.log(x), np.log(_dRdE), 'linear', copy=False, assume_sorted=True, fill_value='extrapolate')
        return np.exp(f(np.log(E)))


class SlabLayer(object):
    """
    SlabLayer(table,material,thickness,units='mils')
    manages slab layer properties and methods
    table - RangeTable or ScreamTable object
    material - string naming material supported by RangeTable
    thickness - numeric thickness
    units - string naming supported units (mils, in, cm, mm, um, g/cm^2)
    """

    def __init__(self, table, material, thickness, units='mils'):
        self.table = table
        self.material = material
        self.thickness = thickness
        self.units = units
        self.density = table.density(material)
        self.cm = self.tocm()
        self.gcm2 = self.cm * self.density
        self._cache = None  # degraded_spectrum cache

    def __repr__(self):
        return 'SlabLayer %s %s %s' % (self.thickness, self.units, self.material)

    def __str__(self):
        return '%s %s %s' % (self.thickness, self.units, self.material)

    def tocm(self, thickness=None, units=None):
        """tocm(thickness=None,units=None)
        thickness is number
        units is string naming any of the units supported by this class
        """
        if thickness is None:
            thickness = self.thickness
        if units is None:
            units = self.units
        if units == 'mils':
            return thickness / 1000 * 2.54
        elif units == 'in':
            return thickness * 2.54
        elif units == 'cm':
            return thickness
        elif units == 'mm':
            return thickness / 10
        elif units == 'um':
            return thickness / 1e4
        elif units == 'g/cm^2':
            return self.tocm(thickness, units) * self.density
        else:
            raise Exception('unknown units %s' % units)

    def Eout(self, Ein, angle=0):
        """
        Eout = Eout(Ein,angle=0)
        Compute exit energy for particle incident with energy Ein in MeV
        Returns exit energy in MeV
        angle in radians from normal
        Returns 0 if absorbed
        expects scalars
        """
        dR = self.gcm2 / np.cos(angle)  # secant adjusted path length through layer
        Rin = self.table.EtoR(self.material, Ein)
        Rout = Rin - dR
        if Rout <= self.table.Rlimits(self.material)[0]:  # below Rlimit, treat as absorbed
            return 0
        else:
            return self.table.RtoE(self.material, Rout)

    def Ein(self, Eout, angle=0):
        """
        Ein = Ein(Eout,angle=0)
        Compute incident energy for particle exiting with energy Eout in MeV
        Returns incident energy in MeV
        angle in radians from normal
        returns inf if above range table Elimits
        """
        dR = self.gcm2 / np.cos(angle)  # secant adjusted path length through layer
        Rout = self.table.EtoR(self.material, Eout)
        Rin = Rout + dR
        if Rin > self.table.Rlimits(self.material)[1]:
            return np.inf
        else:
            return self.table.RtoE(self.material, Rin)

    def dRdE(self, E):
        """
        grad = dRdE(E)
        derivative of range with respect to energy at E
        E in MeV
        grad in g/cm^2/MeV
        """
        return self.table.dRdE(self.material, E)

    def get_cache(self, dEoverE=0.1, epsrel=1e-4):
        """cache = get_cache(dEoverE=0.1,epsrel=1e-4)
        return the cache that stores precomputed transform
        cache['energy'] - cache energy grid (NE,) MeV
        cache['transform'] - transform matrix (NE,NE) dimensionless
        exitflux = dot(transform,influx)
        options:
            dEoverE: relative spacing between cache energy grid points
            epsrel: epsrel argument passed to quad numerical integration (relative precision)
        """
        # build a transform that, for each exiting energy and incident energy, integrates over just the relevant incident angles
        cache = self._cache
        if cache is None:  # prepare cache
            Elimits = self.table.Elimits(self.material)
            NE = np.maximum(100, int(np.ceil(25 * 0.1 / dEoverE * np.log10(Elimits[1] / Elimits[
                0]))))  # number of points in energy integral 100 or enough to have 10% energy spacing
            cache = {'NE': NE}
            E = np.exp(np.linspace(np.log(Elimits[0]), np.log(Elimits[1]), NE))  # energy grid
            cache['energy'] = E
            cache['transform'] = np.zeros((NE, NE))
            # integral is from Emin to infinity dE, instead of 0 to pi/2 dtheta.
            # change of variables leads to
            # Jout = gcm2/2*dRdEout*integral(Emin,inf)dEin Jin/(Rin-Rout)^2
            # Emin is Ein for Eout at normal incidence
            # here's how the change of variables works:
            #  normally, we integrate:
            #  Jout = dRdEout/2*integral(0,pi/2)dtheta*sin(theta)/2*Jin*dRdEin
            #  first, we change x = cos(theta), dx = -sin(theta)dtheta
            #  Jout = dRdEout/2*integral(0,1)dx/2*Jin*dRdEin
            #  Introduce H = inverse of range function, so H(R(E)) = E
            #  new change of variables:
            #  y = Ein(x) = H(Rout + gcm2/x)
            #  x = gcm2/(R(Ein)-Rout) = cos(theta)
            #  dx = -gcm2/(R(Ein)-Rout)^2*dRdEin*dy
            #  Jout = gcm2*dRdEout/2*integral(Emin,inf)Jin/(R(Ein)-Rout)^2dEin
            # the trapezoidal integral replaces Jin with weight that is 0
            # the left and right neighbor, 1 at the grid point, and linear between.
            for iout in range(NE):  # exiting energy index
                Eout = E[iout]
                dRdEout = self.dRdE(Eout)
                Rout = self.table.EtoR(self.material, Eout)
                Emin = self.Ein(Eout)
                for jin in range(iout, NE):  # entering energy index, no lower than exiting energy
                    Ein = E[jin]
                    for side in ['left', 'right']:  # integrate left, right side of trapezoidal rule
                        if side == 'left':
                            if jin == 0:
                                continue  # no neighbor on left side
                            j0 = jin - 1  # index where weight is 0
                        else:  # right
                            if jin == NE - 1:
                                continue  # no neighbor on right side
                            j0 = jin + 1
                        E0 = E[j0]  # energy where weight is 0

                        def func(E):
                            R = self.table.EtoR(self.material, E)
                            dR = R - Rout
                            assert dR >= self.gcm2
                            weight = (E - E0) / (Ein - E0)  # trapezoidal weight: 0 at E0, 1 at Ein
                            # => linearly interpolate spectrum between cache.energy gridpoints
                            return weight / dR ** 2

                        E1 = np.minimum(Ein, E0)
                        E1 = np.maximum(E1, Emin)  # integration cannot start below Emin
                        E2 = np.maximum(Ein, E0)
                        if E2 > E1:
                            (y, abserr) = quad(func, E1, E2, epsrel=epsrel)
                            cache['transform'][iout, jin] += y * self.gcm2 / 2 * dRdEout
            self._cache = cache
        return self._cache

    def degraded_spectrum(self, energy, flux, fast=True, epsrel=1e-4, dEoverE=0.1):
        """
        exitflux = degraded_spectrum(energy,flux,fast=True,epsrel=1e-4,dEoverE=0.1)
        compute degraded spectrum
        energy - array of incident energies
        flux - array of incident fluxes (differential #/cm^2/MeV)
        exitflux - array of transmitted fluxes (differential #/cm^2/MeV)
        options:
            fast - True/False
                if True precomputes transform matrix and uses that on all subsequent calls
                if False performs quad adaptive integral on every call
            epsrel - relative precision of numerical integration, passed to quad
            dEoverE - dEoverE passed to self.get_cache cache energy grid spacing (relative)
        """

        # not Fast:
        # for every exiting energy
        # for every polar angle
        # compute the incident energy
        # and the incident flux at the incident energy
        # then adjust to be the exiting flux
        # exitflux(Eout) = int[0,pi/2] dtheta flux(Ein)*(dR/dEout)/(dR/dEin) sin(theta)/2
        # Ein depends on Eout,theta and thickness+material

        global _HIGH_ENERGY_WARNED  # may modify this is need to warn

        Elimits = self.table.Elimits(self.material)
        if np.all(flux > 0):
            influx = lambda E: np.exp(
                interp1d(np.log(energy), np.log(flux), 'linear', copy=False, assume_sorted=True, bounds_error=False,
                         fill_value=-np.inf)(np.log(E)))
        else:
            influx = lambda E: interp1d(np.log(energy), flux, 'linear', copy=False, assume_sorted=True,
                                        bounds_error=False, fill_value=0)(np.log(E))

        if fast:
            cache = self.get_cache(dEoverE=dEoverE, epsrel=epsrel)
            tmpflux = np.dot(cache['transform'], influx(cache['energy']))
            if np.all(tmpflux > 0):
                exitflux = np.exp(
                    interp1d(np.log(cache['energy']), np.log(tmpflux), 'linear', copy=False, assume_sorted=True,
                             bounds_error=False, fill_value=-np.inf)(np.log(energy)))
            else:
                exitflux = interp1d(np.log(cache['energy']), tmpflux, 'linear', copy=False, assume_sorted=True,
                                    bounds_error=False, fill_value=0)(np.log(energy))
            # copy energies above upper limit, at half intensity
            iHigh = energy > cache['energy'][-1]
            if np.any(iHigh):
                if not _HIGH_ENERGY_WARNED:
                    print('Warning: %g MeV requested. Assuming negligible energy loss for all energies > %g MeV' % (
                    np.min(energy[iHigh]), Elimits[-1]))
                    _HIGH_ENERGY_WARNED = True
                exitflux[iHigh] = flux[iHigh] / 2  # assume highest energies penetrate unaffected
        else:
            exitflux = np.full(flux.shape, np.nan)

            def func(theta, Eout, dRdEout):
                Ein = self.Ein(Eout, theta)
                if (Ein < energy[0]) or (Ein > energy[-1]):
                    return 0
                return influx(Ein) * dRdEout / self.dRdE(Ein) * np.sin(theta) / 2

            Elimits = self.table.Elimits(self.material)
            for (iE, Eout) in enumerate(energy):
                if Eout > Elimits[-1]:
                    if not _HIGH_ENERGY_WARNED:
                        print('Warning: %g MeV requested. Assuming negligible energy loss for all energies > %g MeV' % (
                        Eout, Elimits[-1]))
                        _HIGH_ENERGY_WARNED = True
                    exitflux[iE] = flux[iE] / 2  # assume highest energies penetrate unaffected
                else:
                    dRdEout = self.dRdE(Eout)
                    (y, abserr) = quad(func, 0, np.pi / 2, args=(Eout, dRdEout), epsrel=epsrel)
                    exitflux[iE] = y
        return exitflux


if __name__ == '__main__':
    # perform three demos
    import matplotlib.pyplot as plt

    plt.close('all')

    # load range tables
    stable = ScreamTable()
    ntable = ScreamTable('NISTRangeData_protons.csv')
    etable = ScreamTable('ScreamRangeData_electrons.csv')

    # proton fluence
    pMeV = 10 ** np.linspace(0, 3.5, 51)  # 1 MeV to ~3 GeV
    pFluence = pMeV ** -2  # falling spectrum

    # electron fluence
    eMeV = 10 ** np.linspace(-1, 0.5, 51)  # 0.1 MeV to ~3 MeV
    eFluence = eMeV ** -3  # falling spectrum

    # compare protons SCREAM slow to fast method
    front = SlabLayer(stable, 'SiO2', 15)  # 15 mils of SiO2
    slow = front.degraded_spectrum(pMeV, pFluence, fast=False)
    fast = front.degraded_spectrum(pMeV, pFluence, fast=True)
    plt.figure()
    plt.loglog(pMeV, pFluence / 2, 'k-', pMeV, slow, 'b-', pMeV, fast, 'r--')
    plt.xlabel('MeV')
    plt.ylabel('protons/cm$^2$/MeV')
    plt.title('Comparison of Fast and Slow methods %s' % str(front))
    plt.legend(['Incident/2', 'Slow', 'Fast'])
    plt.grid(True)

    # compare protons SCREAM to NIST (slow method)
    nfront = SlabLayer(ntable, 'SILICON DIOXIDE', 15)  # 15 mils of SiO2
    nslow = nfront.degraded_spectrum(pMeV, pFluence, fast=False)
    plt.figure()
    plt.loglog(pMeV, pFluence / 2, 'k-', pMeV, slow, 'b-', pMeV, nslow, 'r--')
    plt.xlabel('MeV')
    plt.ylabel('protons/cm$^2$/MeV')
    plt.title('Comparison of SCREAM and NIST tables %s' % str(front))
    plt.legend(['Incident/2', 'SCREAM', 'NIST'])
    plt.grid(True)
    plt.text(1e0, 4e-8, 'NIST drops to 0 @ ~3e3 MeV b/c end\n'
                        'of degraded incident spectrum w/in table.\n' + 'SCREAM continues b/c assumes non-degraded\n' + 'spectrum beyond end of table.',
             verticalalignment='bottom', backgroundcolor='w')

    # demo electrons (slow method)
    efront = SlabLayer(etable, 'SiO2', 5)  # 5 mils of SiO2
    eslow = efront.degraded_spectrum(eMeV, eFluence, fast=False)
    plt.figure()
    plt.loglog(eMeV, eFluence / 2, 'k-', eMeV, eslow, 'b-')
    plt.xlabel('MeV')
    plt.ylabel('electrons/cm$^2$/MeV')
    plt.title('Electrons %s' % str(efront))
    plt.legend(['Incident/2', 'Slow Method'])
    plt.grid(True)
