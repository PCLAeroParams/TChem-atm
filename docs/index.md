# **Overview**
TChem-atm computes source terms and Jacobian matrices for chemical systems. It is a performance-portable software toolkit designed for complex kinetic mechanisms. We designed and implemented TChem-atm using [Kokkos](https://github.com/kokkos/kokkos.git).

Software Design:

  * Modern C++.
  * Kokkos programming model for performance portability.
  * CMake build system.
  * Numerical Jacobians and SACADO analytic Jacobians for all models.
  * Coupling to external ODE solvers, e.g., Tines, Sundials (CVODE).

![TChem](figures/TChem_atm.png)

TChem-atm includes a parser for a YAML input file that constructs a kinetic model constant data object containing relevant parameters for the computation of the chemical source terms. It computes the reaction constants and rate of progress for all the reactions listed in the input files. Then, it calculates the net production rate or source terms for all species mentioned in the input files. TChem-atm automatically calculates the Jacobian matrix for the source terms using either finite differences (numerical Jacobian) or automatic differentiation (analytical Jacobian via SACADO). Furthermore, the computation of the source term and associated Jacobian is independent of the time integration solver in TChem-atm. Therefore, TChem-atm provides an interface for time integration (Box model) for the Tines and CVODE libraries. Finally, TChem-atm features a batched interface for evaluating the source term, Jacobian matrix, and time integration.
# **Citations**
* [TChem: A performance portable parallel software
toolkit for complex kinetic mechanisms.](https://www.sciencedirect.com/science/article/pii/S0010465522003472)

```bibtex
@article{tchem:Kim:2022,
  title    = {{TChem: A performance portable parallel software toolkit for complex kinetic mechanisms}},
  journal  = {Computer Physics Communications},
  volume   = {285},
  pages    = {108628},
  year     = {2023},
  issn     = {0010-4655},
  author   = {Kyungjoo Kim and Oscar H. Díaz-Ibarra and Habib N. Najm and Judit Zádor and Cosmin Safta},
  keywords = {TChem, Kokkos, Performance portability, GPU, Flow chemistry}
}
```

* [“Benchmarking TChem for Potential Incorporation into E3SM as a Replacement Chemical Kinetics Solver”](sand_report/QTI_tchemV1.pdf)
```bibtex
@techreport{Diaz-Ibarra:2024:tchem,
  author      = {Diaz-Ibarra, Oscar and Schmidt, Michael J.  and Safta, Cosmin },
  title       = {{Benchmarking TChem for Potential Incorporation into E3SM as a Replacement Chemical Kinetics Solver}},
  institution = {Sandia National Laboratories},
  year        = {2024},
  number      = {SAND2024-01807R}
}
```

# **Installation**
The [installation](installation.md) guide demonstrates how to obtain, build, and install TChem-atm along with its third-party libraries.

# **Theoretical Background**
he TChem-atm approach is briefly described in the [Methodoly section](methodology.md).

# **Input file**
A description of input files is presented in the [Input file section](input.md).

# **Examples**

A list of examples can be found [here](examples.md).

# **Acknowledgements**
TChem-atm has been developed using the following funding sources:

* Sandia Laboratory Directed Research and Development (LDRD) project "Bridging aerosol representations across scales with physics-constrained statistical learning."

* Sandia LDRD project "Benchmarking TChem for Potential Incorporation into E3SM as a Replacement Chemical Kinetics Solver."

* The [EAGLES project](https://climatemodeling.science.energy.gov/projects/enabling-aerosol-cloud-interactions-global-convection-permitting-scales-eagles), which was funded by
the Office of Science's [Biological and Environmental
Research](https://science.osti.gov/ber) Program.

* Exascale Catalytic Chemistry ([ECC](https://www.ecc-project.org/)) Project.
