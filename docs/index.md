# **Overview**
TChem-atm computes source terms and Jacobian matrices for chemical systems. It is performance-portable software toolkit for complex kinetic models and it is designed and implemented using [Kokkos](https://github.com/kokkos/kokkos.git)(a performance portable parallel programming model).

Software Design:

  * Modern C++.
  * Kokkos programming model for performance-portability.
  * CMAKE build system.
  * Numerical Jacobians and SACADO analytic Jacobians for all models.
  * Couples to external ODE solvers. e.g., Tines, Sundials(CVODE).

![TChem](figures/TChem_atm.png)

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
The [Installation](installation.md) guide shows you how to build and install
  TChem-atm on your own machine or on a supported high-performance platform.

# **Theorycal background**
A brief description of TChem-atm approach is presented in [here](methodology.md).

# **Input file**
A description of input files employ is [here](input.md)

# **Examples**

List of examples can be found [here](examples.md).

# **Acknowledgements**
TChem has developed using the following founding sources:

* Sandia Laboratory Directed Research and Development (LDRD) project "Bridging aerosol representations across scales with physics-constrained statistical learning."

* Sandia LDRD project "Benchmarking TChem for Potential Incorporation into E3SM as a Replacement Chemical Kinetics Solver."

* The [EAGLES project](https://climatemodeling.science.energy.gov/projects/enabling-aerosol-cloud-interactions-global-convection-permitting-scales-eagles), which was funded by
the Office of Science's [Biological and Environmental
Research](https://science.osti.gov/ber) Program.

* Exascale Catalytic Chemistry ([ECC](https://www.ecc-project.org/)) Project.
