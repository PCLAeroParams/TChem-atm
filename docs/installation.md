TChem-atm  requires the following third-party libraries:

  * [Tines](https://github.com/sandialabs/Tines.git)
    * [Sacado](https://docs.trilinos.org/dev/packages/sacado/doc/html/index.html)
    * [BLAS/LAPACK](http://www.openblas.net)
    * [YAML](https://github.com/jbeder/yaml-cpp)
    * [Sundials](https://github.com/LLNL/sundials.git)
  * [Gtest](https://github.com/google/googletest.git)
  * [Skywalker](https://github.com/eagles-project/skywalker.git)

These third-party libraries are submodules in TChem-atm. Thus, they can be initialized using:

```bash
git submodule update --init --recursive
```
### **Obtaning TChem-atm**

TChem-atm is open source code available in [GitHub](https://github.com/PCLAeroParams/TChem-atm).

### **Building and installing third-party libraries**

The script [``scripts/tpls_bld.sh``](tpls_bld.sh) clones, builds, and installs the TChem's third-party libraries. To use this script, ones must provide the compiler information using:

```bash
MY_CC=gcc # c++ compiler
MY_CXX=g++ # c++ compiler
MY_FC=gfortran # fortran compiler
```

to build/install with CUDA ``ON`` or ``OFF` using:

```bash
CUDA="OFF" # ON/OFF to compile TChem with NVIDIA-GPUs
```
and the root directory to install/build the third-party libraries;
``ROOT=/path/to/tchem/``

These variables are located in the top section of the ``scripts/tpls_bld.sh``. This ``scripts/tpls_bld.sh`` initializes the submodules and install/build the third-party libries in the ``$ROOT/HOST`` or ``$ROOT/DEVICE`` directory for ``CUDA =OFF`` or ``CUDA=ON``. Because Kokkos can be installed/built for both CUDA ``ON`` and ``OFF``, it will be installed in both director ``$ROOT/HOST`` and  ``$ROOT/DEVICE`` and the other libraries will be installed in  ``$ROOT/HOST``. Hence, one must run this script using CUDA=OFF to install all third party libraries and run it a second time using CUDA=ON, if a NVIDIA GPUs is aviable.

## **Building and install TChem-atm and Tines**
The script [``scripts/build_script.sh``](build_script.sh) builds, and installs the TChem-atm and Tines. Similar to the third-party libraries script, in the ``scripts/build_script.sh`` one must provide the compiler information:
```bash
MY_CC=gcc # c++ compiler
MY_CXX=g++ # c++ compiler
MY_FC=gfortran # fortran compiler
```
to build/install with CUDA ``ON`` or ``OFF`` using:
```bash
CUDA="OFF" # ON/OFF to compile TChem with NVIDIA-GPUs
```
In addition, this script adds the option to turn ON/OFF automatic differenciation using SACADO library using ``SACADO="ON"`` or ``SACADO="OFF"``.
The installation path the host libraries, ``INSTALL_BASE_HOST`` and the ``ROOT=/path/to/tchem``, where TChem-atm and Tines will be installed. If CUDA is ON, libraries will be installed in the ``ROOT/DEVICE`` directory. Otherwise, they will install in the ``ROOOT/HOST`` directory.
