#!/bin/bash

#=======================================================================================
# example script to download, build, install, & test: 
#     openblas, kokkos, gtest, tines, tchem, csplib
#=======================================================================================
# User configuration  -- begin

# example:
# set CUDA="ON" to use GPU, else set to "OFF"
# JFLAG is for makefile compilation on multiple cores, here 4
# if CUDA="ON", make sure nvcc is in your system PATH

MY_CC=gcc-11
MY_CXX=g++-11
MY_FC=gfortran-11
JFLAG="-j 40"
CUDA="OFF" 

# will clone repos under REPO_BASE
# will build under BUILD_BASE
# will install under INSTALL_BASE
# example: as follows under .
ROOT=/Users/odiazib/Documents/p-clap/TChem-atm/external
REPO_BASE=$ROOT/Tines/ext
BUILD_BASE=${PWD}/TPLs/build
INSTALL_BASE=${PWD}/TPLs/install
# User configuration  -- end
#=======================================================================================

#=======================================================================================
# OpenBLAS
# nb. to make sure this openblas gets used when running outside this script
# add in your .bashrc/.bash_profile :
# on mac need this
#    export LIBRARY_PATH="${LIBRARY_PATH}:${OPENBLAS_INSTALL_PATH}/lib"
# on linux need this
#    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${OPENBLAS_INSTALL_PATH}/lib"

build_openblas (){
echo "Building OpenBLAS:"
cd ${OPENBLAS_REPOSITORY_PATH}
make CC=${MY_CC} FC=${MY_FC} HOSTCC=${MY_CC} USE_OPENMP=1
}

install_openblas (){
echo "Installing OpenBLAS:"
cd ${OPENBLAS_REPOSITORY_PATH}
make PREFIX=${OPENBLAS_INSTALL_PATH} install
}

build_install_kokkos(){
echo "Building kokkos:"
mkdir ${KOKKOS_BUILD_PATH}
mkdir ${KOKKOS_INSTALL_PATH}
cd ${KOKKOS_BUILD_PATH}
# kokkos install bin directory is not yet created
cmake \
    -D CMAKE_INSTALL_PREFIX=${KOKKOS_INSTALL_PATH} \
    -D CMAKE_CXX_COMPILER=${KOKKOS_CXX_REPO_COMPILER} \
    -D CMAKE_CXX_FLAGS="-fopenmp -g" \
    -D Kokkos_ENABLE_SERIAL=ON \
    -D Kokkos_ENABLE_OPENMP=ON \
    -D Kokkos_ENABLE_CUDA=${CUDA} \
    -D Kokkos_ENABLE_DEPRECATED_CODE=OFF \
    -D Kokkos_ENABLE_CUDA_CONSTEXPR=${CUDA} \
    -D Kokkos_ENABLE_CUDA_LAMBDA=${CUDA} \
    ${KOKKOS_REPOSITORY_PATH}
make ${JFLAG} install
}
#=======================================================================================

build_install_gtest(){
echo "Building gtest:"
mkdir ${GTEST_BUILD_PATH}
mkdir ${GTEST_INSTALL_PATH}
cd ${GTEST_BUILD_PATH}
cmake \
    -D CMAKE_INSTALL_PREFIX="${GTEST_INSTALL_PATH}" \
    -D CMAKE_CXX_COMPILER="${MY_CXX}"  \
    ${GTEST_REPOSITORY_PATH}
make ${JFLAG} install
}

build_install_yaml(){
echo "Building gtest:"
mkdir ${YAML_BUILD_PATH}
mkdir ${YAML_INSTALL_PATH}
cd ${YAML_BUILD_PATH}
cmake \
    -D CMAKE_INSTALL_PREFIX="${YAML_INSTALL_PATH}" \
    -D CMAKE_CXX_COMPILER="${MY_CXX}"  \
    -D CMAKE_C_COMPILER="${MY_CC}" \
    -D CMAKE_CXX_FLAGS="-g -c" \
    -D CMAKE_EXE_LINKER_FLAGS="" \
    -D CMAKE_BUILD_TYPE=RELEASE \
    ${YAML_REPOSITORY_PATH}
make ${JFLAG} install
}
get_sundials (){
echo "get sundials:"
if [ -d "${SUNDIALS_REPOSITORY_PATH}" ] && [ "$(ls -A ${SUNDIALS_REPOSITORY_PATH})" ]; then
  echo "${SUNDIALS_REPOSITORY_PATH} exists and is not empty ... aborting clone"; return
fi
git clone https://github.com/LLNL/sundials.git ${SUNDIALS_REPOSITORY_PATH}
}

build_install_sundials(){
echo "Building gtest:"
mkdir ${SUNDIALS_BUILD_PATH}
mkdir ${SUNDIALS_INSTALL_PATH}
cd ${SUNDIALS_BUILD_PATH}
cmake \
    -D CMAKE_INSTALL_PREFIX=${SUNDIALS_INSTALL_PATH} \
    -D CMAKE_CXX_COMPILER="${MY_CXX}"  \
    -D CMAKE_C_COMPILER="${MY_CC}" \
    -D CMAKE_CXX_FLAGS="-g" \
    -D CMAKE_C_FLAGS="-g" \
    -D CMAKE_EXE_LINKER_FLAGS="" \
    -D CMAKE_BUILD_TYPE=RELEASE \
    -D ENABLE_CALIPER:BOOL=OFF \
    ${SUNDIALS_REPOSITORY_PATH}
make ${JFLAG} install
}

build_install_skywalker(){
echo "Building skywalker:"
mkdir ${SKYWALKER_BUILD_PATH}
mkdir ${SKYWALKER_INSTALL_PATH}
cd ${SKYWALKER_BUILD_PATH}
cmake \
    -D CMAKE_INSTALL_PREFIX=${SKYWALKER_INSTALL_PATH} \
    -D CMAKE_CXX_COMPILER="${MY_CXX}"  \
    -D CMAKE_C_COMPILER="${MY_CC}" \
    -D SKYWALKER_PRECISION=double \
    -D CMAKE_BUILD_TYPE=RELEASE \
    ${SKYWALKER_REPOSITORY_PATH}
make ${JFLAG} install
}
build_install_tines(){
echo "Building tines:"
mkdir ${TINES_BUILD_PATH}
mkdir ${TINES_INSTALL_PATH}
cd ${TINES_BUILD_PATH}
#    -D TINES_ENABLE_SHARED_BUILD=ON \
cmake \
    -D CMAKE_INSTALL_PREFIX="${TINES_INSTALL_PATH}" \
    -D CMAKE_CXX_COMPILER="${KOKKOS_CXX_COMPILER}" \
    -D CMAKE_CXX_FLAGS="-g" \
    -D CMAKE_C_COMPILER="${MY_CC}" \
    -D CMAKE_EXE_LINKER_FLAGS="-lgfortran" \
    -D TINES_ENABLE_DEBUG=OFF \
    -D TINES_ENABLE_VERBOSE=OFF \
    -D CMAKE_BUILD_TYPE=RELEASE \
    -D TINES_ENABLE_TEST=ON \
    -D TINES_ENABLE_EXAMPLE=ON \
    -D SUNDIALS_INSTALL_PATH="${SUNDIALS_INSTALL_PATH}" \
    -D YAML_INSTALL_PATH="${YAML_INSTALL_PATH}" \
    -D KOKKOS_INSTALL_PATH="${KOKKOS_INSTALL_PATH}" \
    -D GTEST_INSTALL_PATH="${GTEST_INSTALL_PATH}" \
    -D OPENBLAS_INSTALL_PATH="${OPENBLAS_INSTALL_PATH}" \
    ${TINES_REPOSITORY_PATH}/src
make ${JFLAG} install
# if using e.g. yum or apt-get installed openblas, may need this:
# Not needed when using above installed openblas and macport version on osx
#    -D LAPACKE_INSTALL_PATH="${LAPACKE_INSTALL_PATH}" \
}

#=======================================================================================

#=======================================================================================
# main
mkdir -p ${BUILD_BASE}
mkdir -p ${INSTALL_BASE}
OPENBLAS_REPOSITORY_PATH=${REPO_BASE}/OpenBLAS
OPENBLAS_INSTALL_PATH=${INSTALL_BASE}/openblas
#build_openblas
#install_openblas
#
KOKKOS_REPOSITORY_PATH=${REPO_BASE}/kokkos
KOKKOS_BUILD_PATH=${BUILD_BASE}/kokkos
KOKKOS_INSTALL_PATH=${INSTALL_BASE}/kokkos
if [ "${CUDA}" = "ON" ]; then
    KOKKOS_CXX_COMPILER="${KOKKOS_INSTALL_PATH}/bin/nvcc_wrapper"
    KOKKOS_CXX_REPO_COMPILER="${KOKKOS_REPOSITORY_PATH}/bin/nvcc_wrapper"
else
    KOKKOS_CXX_COMPILER="${MY_CXX}"
    KOKKOS_CXX_REPO_COMPILER="${MY_CXX}"
fi
#build_install_kokkos
#
GTEST_REPOSITORY_PATH=${REPO_BASE}/gtest
GTEST_BUILD_PATH=${BUILD_BASE}/gtest
GTEST_INSTALL_PATH=${INSTALL_BASE}/gtest
#build_install_gtest
#
YAML_REPOSITORY_PATH=${REPO_BASE}/yaml
YAML_BUILD_PATH=${BUILD_BASE}/yaml
YAML_INSTALL_PATH=${INSTALL_BASE}/yaml
#build_install_yaml

SUNDIALS_REPOSITORY_PATH=$ROOT/sundials
SUNDIALS_BUILD_PATH=${BUILD_BASE}/sundials
SUNDIALS_INSTALL_PATH=${INSTALL_BASE}/sundials
#get_sundials
#build_install_sundials

SKYWALKER_REPOSITORY_PATH=$ROOT/Skywalker
SKYWALKER_BUILD_PATH=${BUILD_BASE}/skywalker
SKYWALKER_INSTALL_PATH=${INSTALL_BASE}/skywalker
#build_install_skywalker

TINES_REPOSITORY_PATH=$ROOT/Tines
TINES_BUILD_PATH=${BUILD_BASE}/tines
TINES_INSTALL_PATH=${INSTALL_BASE}/tines
build_install_tines

exit
