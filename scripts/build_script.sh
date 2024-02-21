#!/bin/bash

#=======================================================================================
# User configuration  -- begin

MY_CC=gcc-11
MY_CXX=g++-11
MY_FC=gfortran-11
JFLAG="-j 40"
CUDA="OFF"

INSTALL_BASE_HOST=/Users/odiazib/Documents/p-clap/build_tchem/TPLs/install
SKYWALKER_INSTALL_PATH=${INSTALL_BASE_HOST}/skywalker
OPENBLAS_INSTALL_PATH=${INSTALL_BASE_HOST}/openblas
GTEST_INSTALL_PATH=${INSTALL_BASE_HOST}/gtest
YAML_INSTALL_PATH=${INSTALL_BASE_HOST}/yaml
SUNDIALS_INSTALL_PATH=${INSTALL_BASE_HOST}/sundials

if [ "${CUDA}" = "ON" ]; then
REPO_BASE=${ROOT}/CSPlib/ext/TChem/ext/Tines/ext
BUILD_BASE=${PWD}/DEVICE/RELEASE/build
INSTALL_BASE=${PWD}/DEVICE/RELEASE/install
INSTALL_BASE_DEVICE=${PWD}/DEVICE/install
KOKKOS_REPOSITORY_PATH=${REPO_BASE}/kokkos
KOKKOS_INSTALL_PATH=${INSTALL_BASE_DEVICE}/kokkos
KOKKOS_CXX_COMPILER="${KOKKOS_INSTALL_PATH}/bin/nvcc_wrapper"
KOKKOS_CXX_REPO_COMPILER="${KOKKOS_REPOSITORY_PATH}/bin/nvcc_wrapper"
else
BUILD_BASE=${PWD}/HOST/RELEASE/build
INSTALL_BASE=${PWD}/HOST/RELEASE/install
KOKKOS_INSTALL_PATH=${INSTALL_BASE_HOST}/kokkos
TINES_INSTALL_PATH=${INSTALL_BASE_HOST}/tines
KOKKOS_CXX_COMPILER="${MY_CXX}"
KOKKOS_CXX_REPO_COMPILER="${MY_CXX}"
fi

# User configuration  -- end

#=======================================================================================
build_install_tchem_atm(){
echo "Building TChem_atm:"
mkdir ${TCHEM_BUILD_PATH}
mkdir ${TCHEM_INSTALL_PATH}
cd ${TCHEM_BUILD_PATH}
#    -D TCHEM_ENABLE_SHARED_BUILD=ON \
cmake \
    -D CMAKE_INSTALL_PREFIX="${TCHEM_INSTALL_PATH}" \
    -D CMAKE_CXX_COMPILER="${KOKKOS_CXX_COMPILER}" \
    -D CMAKE_C_COMPILER="${MY_CC}" \
    -D CMAKE_EXE_LINKER_FLAGS="-lgfortran" \
    -D CMAKE_CXX_FLAGS="-g" \
    -D TCHEM_ATM_ENABLE_VERBOSE=OFF \
    -D CMAKE_BUILD_TYPE=RELEASE \
    -D TCHEM_ATM_ENABLE_TEST=ON \
    -D TCHEM_ATM_ENABLE_EXAMPLE=ON \
    -D TCHEM_ATM_ENABLE_SACADO_JACOBIAN_ATMOSPHERIC_CHEMISTRY=OFF\
    -D KOKKOS_INSTALL_PATH="${KOKKOS_INSTALL_PATH}" \
    -D TINES_INSTALL_PATH="${TINES_INSTALL_PATH}" \
    -D TCHEM_ATM_ENABLE_SKYWALKER=ON \
    -D SKYWALKER_INSTALL_PATH=${SKYWALKER_INSTALL_PATH} \
    -D GTEST_INSTALL_PATH="${GTEST_INSTALL_PATH}" \
    ${TCHEM_REPOSITORY_PATH}/src
make ${JFLAG} install
}
#=======================================================================================
#=======================================================================================
# main


mkdir -p ${BUILD_BASE}
mkdir -p ${INSTALL_BASE}

TCHEM_REPOSITORY_PATH=/Users/odiazib/Documents/p-clap/TChem-atm/
TCHEM_BUILD_PATH=${BUILD_BASE}/tchem_atm
TCHEM_INSTALL_PATH=${INSTALL_BASE}/tchem_atm
build_install_tchem_atm

exit
