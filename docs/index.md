# Overview
TChem-atm computes source terms and Jacobian matrices for chemical systems. It is performance-portable software toolkit for complex kinetic models and it is designed and implemented using [Kokkos](https://github.com/kokkos/kokkos.git)(a performance portable parallel programming model).

![TChem](figures/TChem_atm.png)

Software Design:

  * Modern C++.
  * Kokkos programming model for performance-portability.
  * CMAKE build system.
  * Numerical Jacobians and SACADO analytic Jacobians for all models.
  * Couples to external ODE solvers. e.g., Tines, Sundials(CVODE).

# Installation
The [Installation](installation.md) guide shows you how to build and install
  TChem-atm on your own machine or on a supported high-performance platform.

# Theorycal background
A brief description of TChem-atm approach is presented in [here](methodology.md).

# Input file
A description of input files employ is [here](input.md)

## Tests

## Examples

## Implimentation

## Acknowledgements
TChem has developed using the following founding sources:
* LDRD sandia
* The [EAGLES project](https://climatemodeling.science.energy.gov/projects/enabling-aerosol-cloud-interactions-global-convection-permitting-scales-eagles), which was funded by
the Office of Science's [Biological and Environmental
Research](https://science.osti.gov/ber) Program.
* ECC

The source code is available on
[GitHub](https://github.com/PCLAeroParams/TChem-atm).
