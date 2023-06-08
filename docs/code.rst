EQFLUX
========================================
This code is taken almost line for line from the JPL GaAs Solar Cell Radiation Handbook.  The handbook publishes the FORTRAN code for EQGAFLUX.  EQGAFLUX takes RDC data, glass range tables for protons and electrons, GaAs stopping table, and an input space environment particle spectrum to determine the 1 MeV equivalent electron fluence.  This module takes that original FORTRAN code and replicates it line for line in python code.

.. automodule:: EQFLUX
    :members:


degradation_equations
========================================
A set of what is know as the degradation equation as described in the JPL Solar Cell and GaAs Solar Cell Handbook.  The equation is used to fit, interpolate, and extrapolate solar cell radiation degrdation data to generate relative damage coefficients as well as displacement damage dose curves.  These equations have been empirically determined to best approximate radiation degradation data, but they do not always fit well.  This is an active area to determine the best way to relate solar cell radiation degradation to fundamental physics that govern solar cell degradatio due to radiation

.. automodule:: degradation_equations
    :members:

electron_rdc
========================================
Relative damage coefficients (RDC) are used to relate the relative damage of one particle energy to another.  The JPL radiation handbooks relate all particle energies and types to a 1MeV electron equivalent fluence, hence the name JPL Equivalent Fluence of JPL EQFLUX method.  According the JPL Handbooks, 1 MeV electron was choose because it is readily available for ground testing and the abundance of 1 MeV electron particles in the space environment.

Note:  This module and the proton module can be combined to one rdc module.  The modules are seperate to reflect the procedure used in the JPL radiation hand books.

.. automodule:: electron_rdc
    :members:

proton_rdc
========================================
Relative damage coefficients (RDC) are used to relate the relative damage of one particle energy to another.  The JPL radiation handbooks relate all particle energies and types to a 1MeV electron equivalent fluence, hence the name JPL Equivalent Fluence of JPL EQFLUX method.  According the JPL Handbooks, 1 MeV electron was choose because it is readily available for ground testing and the abundance of 1 MeV electron particles in the space environment. For protons all protons are related to 10 MeV protons.  The 10 MeV protons are then related to 1 MeV electrons.

Note:  This module and the proton module can be combined to one rdc module.  The modules are seperate to reflect the procedure used in the JPL radiation hand books

Note: Turns out since JPL EQFLUX was implemented and through DDD 10 MeV protons don't actually do much damage......sooooo hmmmmm would it be better to relate all protons to a lower proton energy.

.. automodule:: proton_rdc
    :members:

relativeMeV_fluence
========================================
This module takes the relative damage coefficients of a solar cell type and the space environment particle spectrum to arrive and total 1 MeV fluence.  Using the 1 MeV electrion radiation data used to generate the RDC curve, one can predict the expected remaining factor of the parameter of interest in a solar cell

.. automodule:: relativeMeV_fluence
    :members:

ddd
========================================
Displacement damage dose (DDD) tools to convert solar cell radiation data to displacement damage dose using the non-ionizing energy loss (NIEL) for a solar cell material of interest.  DDD uses the NIEL curve to determine the displacement damage dose from each particle type and energy as opposed to the RDC method in the JPL EQFLUX method.

.. automodule:: ddd_degradation_functions
    :members:

Example Class
========================================
Several classes built with the library to simplify use. 

.. automodule:: eqflux_aero
    :members:
