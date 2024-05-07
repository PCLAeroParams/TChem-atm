## **UCI mechanism - E3SM version 3 chemistry**

E3SM version 3 uses the UCI mechanism to model gas chemistry and gas - aerosol interations. In E3SM version 3 code, two separed reaction mechanisms are employed one for the Troposphere and Stratosphere. As part the Sandia LDRD project "Benchmarking TChem for Potential Incorporation into E3SM as a Replacement Chemical Kinetics Solver", we created input files and implemented the reaction types that are needed to reproduce this chemistry mechanism in TChem-atm.

E3SM employes the Community Atmospheric Model Pre-Processor (CAMPP) to generate a set of Fortran files to represent and solve the chemical kinetic model, which includes reaction coefficients, the hand side of volumetric mixing ratio (vmr, denoted $\eta$), and ODE solvers (Implicit/Explicit solver). Ultimately, the goal of CAMPP is to compute an updated vmr ($\eta_{t+\Delta t}$) for the troposphere and stratosphere.

For aditional information on this LDRD project, we recommend to read the [SAND2024-01807R report](/../sand_report/QTI_tchemV1.pdf). As we described in the SAND2024-01807R report, we employed TChem-atm to time advance one time step and compare this outputs with E3SM output in six different locations; see table below. We compute the relative root mean square error (RRMSE) of the difference between the E3SM and TChem-atm outputs. Using this approach, we were able to verify the input files and implementation of TChem-atm.

| Abbreviation | Latitude | Longitude | Location |
|--------------|----------|-----------|----------|
| LA | 34.0549 | -118.2426 | Los Angeles, CA, USA |
| BRW | 71.323 | -156.6114 | North Slope, AK, USA |
| MHD | 53.326 | -9.899 | Halfmace, County Galway, Ireland |
| PSA | -64.7742 | -64.0527 | Palmer Station, Antarctica |
| RPB | 13.165 | -59.432 | Ragged Point, Barbados |
| SYO | -69.0125 | 39.59 | Showa Station, Antarctica |
| ZEP | 78.9067 | 11.8883 | Zeppelin mountain, Ny-\AA lesund, Norway |


### **Troposphere mechanism**
The bash scripts and input files for the UCI chemistry in the troposhere are located at ``src/examples/runs/uci_col``. In this directory, we split the input file of TChem in two yaml files. One with initial conditions, and the other file with the list of reactions, constant species and active species. The chemistry file of the UCI mechanism is saved at ``src/examples/runs/uci_col/uci_v2_test3.yaml``. The initial conditions for one UCI test is save at ``src/examples/runs/uci_col/input_conditions_col.yaml`, which corresponds a set of cells for one location in a E3SM simulation.

Mechanism details:

* Number of species : 82
* Number of invariants or constant species: 29
* Number of reaction : 106
* Number of Arrhenius type reactions : 71
* Number of JPL-Troe type reactions : 10
* Number of Ratio JPL-Arrhenius type reactions : 3
* Number of Photolysis rates : 22

A run script for the UCI mechanism is presented next:

```bash
exec=$TCHEM_INSTALL_PATH/examples/TChem_AtmosphericChemistryE3SM.implicit_euler.x
input=$TCHEM_INSTALL_PATH/examples/runs/atmospheric_chemistry/uci_col/uci_v2_test3.yaml
inputFile=$TCHEM_INSTALL_PATH/examples/runs/atmospheric_chemistry/uci_col/input_conditions_multi_col.yaml
run_this="$exec --chemfile=$input \
          --inputfile=$inputFile \
          --outputfile=full_gas.dat \
          --time-iterations-per-interval=100 \
          --tol-time=1e-6 \
          --dtmin=1800 \
          --dtmax=1800 \
          --tend=1800 \
          --atol-newton=1e-18 \
          --rtol-newton=1e-8 \
          --max-newton-iterations=20 \
          --max-time-iterations=20000"

echo $run_this
eval $run_this
```
Here, TChem-atm is using an implicit euler solver. The ``inputfile`` flag allows to pass the initial condition files as an independent file of the ``chemfile``. Furthermore, TChem-atm also offers executables for Tines-TrBDF2 solver, ``TChem_AtmosphericChemistryE3SM.x`` and for CVODE, ``TChem_AtmosphericChemistryE3SM_CVODE.x``; Note the CVODE-TChem-atm executable works only in CPUs.

![E3SM vs TChem-atm Troposhere](figures/net_production_rate_worst.png)
Parity plot for rate of progress in the troposphere. E3SM outputs are produced by CAMPP’s generated code. There are 104 reactions, and we only display the 10 with the larges RRMSE for the locations given in previous table. RRMSE per location is presented in each plot. The net production rates correspond to the RHS of the equations solved by E3SM and TChem.


### **Stratosphere mechanism**

The chemistry file of the UCI mechanism of the stratosphere is saved at ``src/examples/runs/uci_col/uci_explicit_mech.yaml``. The initial conditions for one UCI test is save at ``src/examples/runs/uci_col/input_conditions_explicit_part_multi_col.yaml`, which corresponds to one set of cells for one location in a E3SM simulation.

Mechanism details:

* Number of species : 7
* Number of invariants or constant species invariants : 3
* Number of reaction : 5
* Number of Arrhenius type reactions : 3
* Number of JPL-Troe type reactions : 2

E3SM uses an explicit euler solver to time integrate gas chemistry in the stratoshere. In TChem-atm, we also implimented an explicit euler solver.

```bash
exec=$TCHEM_INSTALL_PATH/examples/TChem_AtmosphericChemistryE3SM.explicit_euler.x
input=$TCHEM_INSTALL_PATH/examples/runs/atmospheric_chemistry/uci_col/uci_explicit_mech.yaml
inputFile=$TCHEM_INSTALL_PATH/examples/runs/atmospheric_chemistry/uci_col/input_conditions_explicit_part_multi_col.yaml
run_this="$exec --chemfile=$input \
          --inputfile=$inputFile \
          --outputfile=full_gas_stratosphere.dat \
          --time-iterations-per-interval=100 \
          --tol-time=1e-6 \
          --dtmin=1800 \
          --dtmax=1800 \
          --tend=1800 \
          --atol-newton=1e-18 \
          --rtol-newton=1e-8 \
          --max-newton-iterations=20 \
          --max-time-iterations=20000"

echo $run_this
eval $run_this
```

![E3SM vs TChem-atm Stratoshere](figures/net_production_rates_stratosphere.png)

Parity plot of net production rates(RHS) in the stratosphere. RRMSE per location is presented in each plot. The net production rates correspond to the RHS of the equations solved by E3SM and TChem.
