# Cronos Celerity

Port of magnetohydrodynamics simulation software `Cronos` using `Celerity`.

## Usage

### Prerequisites

- Celerity
- C++20 capable backend compiler
    - DPCPP works with latest mainline release of Celerity
    - AdaptiveCPP only works with the instruction graph branch (`instruction-graph-wip`) of Celerity
- cmake

### Compiling with DPCPP

1. `cmake -DCMAKE_CXX_COMPILER=/path/to/dpcpp/lib/bin/clang++ -GNinja -Bbuild -DCMAKE_PREFIX_PATH="/path/to/celerity" -DCMAKE_BUILD_TYPE=Release -DCELERITY_DPCPP_TARGETS="nvptx64-nvidia-cuda"`
2. `ninja -C build`

Note: The `lib` folder of dpcpp must be in `LD_LIBRARY_PATH`:

```console
$ export LD_LIBRARY_PATH=/path/to/dpcpp/lib/lib:$LD_LIBRARY_PATH
```

### Compiling with AdaptiveCPP

1. `cmake -GNinja -Bbuild -DCMAKE_PREFIX_PATH="/path/to/celerity" -DCMAKE_BUILD_TYPE=Release -DHIPSYCL_TARGETS="cuda:sm_75" -DUSE_OpenMP=OFF`
2. `ninja -C build`

### Run

Execute `./cronos-celerity` in the build directory, with the full or relative path to
a configuration file as command line argument. Examples can be found in the
folder `configuration/`.

To run the CPU implementation, use `./cronos-amr`, which in the default configuration uses the same runge-kutta solver as the celerity version. The CPU version additionally also supports an adaptive mesh refinement implementation which can be enabled through the config with `enable_amr`.

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

To run the tests, execute the `./tests` binary.

#### Benchmarking

There is a benchmarking binary `./cronos-benchmark` which takes the implementation to be benchmarked (one of `cpu|sycl|celerity`) and the problem size. The problem size is a number which expresses how many times the default size `4x1x1` should be doubled.
