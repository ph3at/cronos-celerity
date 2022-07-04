# Cronos AMR

Reimplementation of magnetohydrodynamics simulation software `Cronos` using
`Adeptive Mesh Refinement (AMR)`.

## Usage

### Compile

1. `cmake -B build -DCMAKE_BUILD_TYPE=Release`
2. `cd build`
3. `make`

### Run

Execute `./cronos-amr` in the build directory, with the full or relative path to
a configuration file as command line argument. Examples can be found in the
folder `configuration/`.

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
need to be registered in `src/user-interface/config-parser.cpp` (line 13).

### Test

To run the tests, simply use `ctest`.
