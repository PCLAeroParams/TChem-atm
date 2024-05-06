## **Carbon Bound 5**
We adapted the carbon bound 5 reaction mechanism from the [CAMP chemistry code](https://github.com/open-atmos/camp/tree/main/mechanisms/cb05cl_ae5).

Mechanism details:

* Number of Species : 67
* Number of Reactions : 187
* Arrhenius type reactions : 168
* CMAQ_H2O2 type reactions : 3
* Troe type reactions : 16
* Sources e.g., EMISSION : 14

The mechanism in CAMP solver has 26 photolysis reactions, which we replace with Arrhenius type reaction with energy = 0  and temperature coefficient = 1.

Scripts to run and plots the outputs of this examples are at : ``src/examples/runs/atmospheric_chemistry/CB05CL_AE5``. The bash script to run this mechanism is :


```bash
exec=$TCHEM_INSTALL_PATH/examples/TChem_AtmosphericChemistry.x

run_this="$exec --chemfile=config_full_gas.yaml \
          --outputfile=full_gas.dat \
          --time-iterations-per-interval=10 \
          --tol-time=1e-10 \
          --dtmin=1e-20 \
          --dtmax=10 \
          --tend=600\
          --atol-newton=1e-18 \
          --rtol-newton=1e-8 \
          --max-newton-iterations=20 \
          --max-time-iterations=20000"

echo $run_this
eval $run_this
```

Here, the ``TChem_AtmosphericChemistry.x`` executable is a box model that time integrates a list of species using the mechanism file from input ``chemfile``, which in this case is the ``config_full_gas.yaml ``. In this box model, only chemical reaction are considered. The outputs of this executable are saved in ``outputfile=full_gas.dat``. In the example directory, we have provided the jupyter-notebook ``PlotFullGas``, which plots the time profiles of each species from the TChem-atm and CAMP outputs. Note, that the CAMP output were previous computed and saved in the TChem-atm repo.

The ``TChem_AtmosphericChemistry.x`` executable allows to change time integration parameters of the [TrBDF2 solver](https://github.com/sandialabs/Tines?tab=readme-ov-file#timeintegration). First, the newton solver parameters are the
the absolute (``atol-newton``), relative tolerance (``rtol-newton``), and the max number of newton iterations (``max-newton-iterations``).  Second, the time step size is controlled using the ``tol-time`` parameter and the maximun(``dtmax``) and minimum (``dtmin``) time step. Third, the parameters ``tend`` or ``max-time-iterations`` will end the simulation. Finally, one can get aditional help from the ``TChem_AtmosphericChemistry.x`` executable using ``TChem_AtmosphericChemistry.x --help``.


![TChem-atm vs CAMP](figures/carbonB5TChemvsCAMP.png)
Comparing of TChem-atm and CAMP outputs using the carbon bound 5 mechanism.
