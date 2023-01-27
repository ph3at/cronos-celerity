# Cronos SYCL

Port of magnetohydrodynamics simulation software `Cronos` using `SYCL`.

## Usage

### Prerequisites

- hipSYCL with C++20 capable backend compiler
- cmake

### Compiling

1. `cmake -Bbuild -DCMAKE_PREFIX_PATH="/path/to/hipSYCL" -DHIPSYCL_TARGETS="cuda:sm_52" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=clang++`
2. `cd build`
3. `make -j`

### Run

Execute `./cronos-sycl` in the build directory, with the full or relative path to
a configuration file as command line argument. Examples can be found in the
folder `configuration/`.

To run the CPU implementation, use `./cronos-amr`, which in the default configuration uses the same runge-kutta solver as the sycl version. The CPU version additionally also supports an adaptive mesh refinement implementation which can be enabled through the config with `enable_amr`.

### Configure

Problem parameters are given by configuration files in the TOML format. To change
the configuration, simply write your own file or (copy and) edit the examples.
The configuration file `shock-tube-integration.toml` needs to be kept unchanged,
as it contains the parameters for the integration test.

`problem_type` is the only required field in configuration files. It can only use
existing problem types. New ones can be created by sub-classing the `Problem`
class, and implementing `initialiseGrid`, `applyBoundary`, and `applySource`,
as well as the `copy` and `parser constructors`. For reference, take a look at
`src/configuration/shock-tube.h` and `-.cpp`. Furthermore, new problem types
need to be registered in `src/user-interface/config-parser.cpp`.

### Test

To run the tests, simply use `ctest`.

#### Benchmarking

The tests also include a benchmarking test for CPU and GPU which automatically scales the problem size by a factor of two. By default the benchmark does 5 doubling steps and repeats each run 3 times.
