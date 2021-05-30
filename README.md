# CDM

The comprehensive Dynamics motion (CDM) is a library built on top of COMA library to run algorithms along N-CMTM tool.

## Installation

The V1 is still in WIP and not all algorithms are working properly.
The V0 is an executable that has working algorithms, tests and benchmarks.

### Dependencies

To compile you need the following tools:

 * [Git]()
 * [CMake]() >= 3.8.2
 * [doxygen]()
 * [g++]() >= 8.3 (for C++17 support)
 * [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.2
 * [coma](https://github.com/vsamy/coma)
 * [RBDyn](https://github.com/jrl-umi3218/RBDyn)

For v0:

 * [CppAD](https://github.com/coin-or/CppAD)
 * [benchmark](https://github.com/google/benchmark)

### Building

```sh
git clone --recursive https://github.com/vsamy/coma
cd cdm
mkdir build
cd build
cmake [options] ..
make -j8 install
```

### Building V0

```sh
git clone --recursive https://github.com/vsamy/coma
cd cdm
mkdir build
cd build
cmake -DBUILD_V0=ON ..
make -j8
```