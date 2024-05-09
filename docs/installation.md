TChem-atm requires the following third-party libraries:

  * [Tines](https://github.com/sandialabs/Tines.git)
    * [Sacado](https://docs.trilinos.org/dev/packages/sacado/doc/html/index.html)
    * [BLAS/LAPACK](http://www.openblas.net)
    * [YAML](https://github.com/jbeder/yaml-cpp)
    * [Sundials](https://github.com/LLNL/sundials.git)
  * [Gtest](https://github.com/google/googletest.git)
  * [Skywalker](https://github.com/eagles-project/skywalker.git)

These third-party libraries are submodules in TChem-atm. Thus, we can initialize them using the following command:

```bash
git submodule update --init --recursive
```
### **Obtaning TChem-atm**
<!--We will need to update the gihub link afther the library is realeased.-->
TChem-atm is open-source code available on [GitHub](https://github.com/PCLAeroParams/TChem-atm).

### **Building and installing third-party libraries**

The script ``scripts/tpls_bld.sh`` clones, builds, and installs TChem-atm's third-party libraries. To use this script, one must provide compiler information as follows:

```bash
MY_CC=gcc # C++ compiler
MY_CXX=g++ # C++ compiler
MY_FC=gfortran # Fortran compiler
```
To build/install with CUDA ``ON`` or ``OFF``, use:

``
CUDA # Set to ON/OFF to compile TChem-atm with or without NVIDIA-GPUs support.
``
Additionally, specify the root directory for installing/building the third-party libraries:

``ROOT=/path/to/tchem-atm/``.

These variables are located in the top section of the ``scripts/tpls_bld.sh`` script.

This script initializes the submodules and installs/builds the third-party libraries in the ``$ROOT/HOST`` or ``$ROOT/DEVICE`` directory, depending on whether ``CUDA=OFF`` or ``CUDA=ON``. Since Kokkos can be installed/built with both CUDA ``ON`` and ``OFF``, this script installs it in both ``$ROOT/HOST`` and ``$ROOT/DEVICE``, while the other libraries are installed in ``$ROOT/HOST``. Therefore, one must run this script using ``CUDA=OFF`` to install all third-party libraries and run it a second time with ``CUDA=ON`` if an NVIDIA GPU is available.

## **Building and install TChem-atm and Tines**
The script [``scripts/build_script.sh``](build_script.sh) builds and installs the TChem-atm and Tines. Similar to the third-party libraries script, in the ``scripts/build_script.sh`` one must provide the compiler information:
```bash
MY_CC=gcc # c++ compiler
MY_CXX=g++ # c++ compiler
MY_FC=gfortran # fortran compiler
```
to build/install with CUDA ``ON`` or ``OFF`` using:
```bash
CUDA="OFF" # ON/OFF to compile TChem with NVIDIA-GPUs
```
In addition, this script adds the option to turn ON/OFF automatic differentiation using the SACADO library using ``SACADO="ON"`` or ``SACADO="OFF"``.
The installation-path the third-party libraries, ``INSTALL_BASE_HOST`` and the ``ROOT=/path/to/tchem``, where TChem-atm and Tines are installed. If CUDA is ON, libraries will be installed in the ``ROOT/DEVICE`` directory. Otherwise, they will install in the ``ROOT/HOST`` directory.
