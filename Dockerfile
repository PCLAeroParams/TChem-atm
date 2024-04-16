FROM ubuntu:22.04

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        autoconf \
        clang-format \
        cmake \
        gcc \
        g++ \
        gfortran \
        git \
        lcov \
        make \
        libopenblas-dev \
        pkg-config \
        ca-certificates

COPY . /tchem_dir/

RUN cmake -S /tchem_dir/external/Tines/ext/kokkos -B /build/kokkos_build \
          -DCMAKE_INSTALL_PREFIX="/install/kokkos_install" \
          -DCMAKE_BUILD_TYPE=RELEASE \
          -DCMAKE_CXX_COMPILER=g++ \
          -DCMAKE_CXX_FLAGS="-fopenmp -g" \
          -DKokkos_ENABLE_SERIAL=ON \
          -DKokkos_ENABLE_OPENMP=ON \
          -DKokkos_ENABLE_CUDA=OFF \
          -DKokkos_ENABLE_CUDA_CONSTEXPR=OFF \
          -DKokkos_ENABLE_CUDA_LAMBDA=OFF
     
WORKDIR /build/kokkos_build/
RUN make -j \
    && make install

RUN cmake -S /tchem_dir/external/Tines/ext/gtest -B /build/gtest_build \
          -DCMAKE_INSTALL_PREFIX="/install/gtest_install" \
          -DCMAKE_CXX_COMPILER=g++
WORKDIR /build/gtest_build
RUN make -j \
    && make install

RUN cmake -S /tchem_dir/external/Tines/ext/yaml -B /build//yaml_build \
          -DCMAKE_INSTALL_PREFIX="/install/yaml_install" \
          -DCMAKE_CXX_COMPILER=g++ \
          -DCMAKE_C_COMPILER=gcc \
          -DCMAKE_CXX_FLAGS="-g -c" \
          -DCMAKE_EXE_LINKER_FLAGS="" \
          -DCMAKE_BUILD_TYPE=RELEASE
WORKDIR /build/yaml_build
RUN make -j \ 
    && make install

RUN cmake -S /tchem_dir/external/Sundials -B /build/sundials_build \
          -DCMAKE_INSTALL_PREFIX="/install/sundials_install" \
          -DCMAKE_CXX_COMPILER=g++ \
          -DCMAKE_C_COMPILER=gcc \
          -DCMAKE_CXX_FLAGS="-g" \
          -DCMAKE_C_FLAGS="-g" \
          -DCMAKE_EXE_LINKER_FLAGS="" \
          -DENABLE_CALIPER:BOOL=OFF \
          -DCMAKE_BUILD_TYPE=RELEASE
WORKDIR /build/sundials_build
RUN make -j \
    && make install

RUN cmake -S /tchem_dir/external/Skywalker -B /build/skywalker_build \
          -DCMAKE_INSTALL_PREFIX="/install/skywalker_install" \
          -DCMAKE_CXX_COMPILER=g++ \
          -DCMAKE_C_COMPILER=gcc \
          -DSKYWALKER_PRECISION=double \
          -DCMAKE_BUILD_TYPE=RELEASE
WORKDIR /build/skywalker_build
RUN make -j \
    && make install

RUN cmake -S /tchem_dir/external/Tines/src -B /build/tines_build \
          -DCMAKE_INSTALL_PREFIX="/install/tines_install" \
          -DCMAKE_CXX_COMPILER=g++ \
          -DCMAKE_CXX_FLAGS="-g" \
          -DCMAKE_C_COMPILER=gcc \
          -DCMAKE_EXE_LINKER_FLAGS="-lgfortran" \
          -DTINES_ENABLE_DEBUG=OFF \
          -DTINES_ENABLE_VERBOSE=OFF \
          -DCMAKE_BUILD_TYPE=RELEASE \
          -DTINES_ENABLE_TEST=OFF \
          -DTINES_ENABLE_EXAMPLE=OFF \
          -DSUNDIALS_INSTALL_PATH=/install/sundials_install \
          -DYAML_INSTALL_PATH=/install/yaml_install \
          -DKOKKOS_INSTALL_PATH=/install/kokkos_install \
          -DOPENBLAS_INSTALL_PATH=`/usr/lib64` \
          -DGTEST_INSTALL_PATH=/install/gtest_install
WORKDIR /build/tines_build
RUN make -j \
    && make install

RUN cmake -S /tchem_dir/src -B /tchem_build \
          -DCMAKE_INSTALL_PREFIX=/tchem_install \
          -DCMAKE_CXX_COMPILER=g++ \
          -DCMAKE_CXX_FLAGS="-g" \
          -DCMAKE_C_COMPILER=gcc \
          -DCMAKE_EXE_LINKER_FLAGS="-lgfortran" \
          -DCMAKE_BUILD_TYPE=RELEASE \
          -DTCHEM_ATM_ENABLE_VERBOSE=OFF \
          -DTCHEM_ATM_ENABLE_TEST=ON \
          -DTCHEM_ATM_ENABLE_EXAMPLE=ON \
          -DTCHEM_ATM_ENABLE_SACADO_JACOBIAN_ATMOSPHERIC_CHEMISTRY=ON\
          -DKOKKOS_INSTALL_PATH=/install/kokkos_install \
          -DTINES_INSTALL_PATH=/install/tines_install \
          -DTCHEM_ATM_ENABLE_SKYWALKER=ON \
          -DTCHEM_ATM_ENABLE_REAL_TYPE="double" \
          -DSKYWALKER_INSTALL_PATH=/install/skywalker_install \
          -DGTEST_INSTALL_PATH=/install/gtest_install
WORKDIR /tchem_build
RUN make -j \
    && make install
