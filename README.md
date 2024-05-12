# arora-ge-ntru

A C++ implementation of the Arora-Ge algorithm for breaking NTRU with multiple keys.

There is also a SageMath implementation in the `sagemath` directory.

This software was developed to supplement the article
"An algebraic algorithm for breaking NTRU with multiple keys" by
Shi Bai, Hansraj Jangir, Tran Ngo, and William Youmans.

# Requirements
The [FLINT](https://github.com/flintlib/flint) library is required.
This program can see some benefit if FLINT is compiled with BLAS enabled, 
for example with [OpenBLAS](https://www.openblas.net/).
[CMake](https://cmake.org/download/) is also required.


# Installation
Install with:
```
mkdir build && cd build
cmake ..
make
```

If you installed FLINT (and its dependencies) in a location other than `/usr/local`
you may need to specify the directory to cmake, for example:
```
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local
```

# Quick start
Run `quickstart.sh` in the top directory to run all steps of the algorithm with
small default parameters, which can be changed on the command line.

# Usage
The main binary is located at `build/apps/bin/arora_ge_ntru`.
Run it with no arguments for an explanation of how to use it.
```
$ ./arora-ge-ntru
Usage: arora-ge-ntru [--help] [--version] [--coeffs VAR] [--seed VAR] [--ring VAR] [--verbose] n q {all,keygen,recover,system,verify}

Arora-Ge algorithm for NTRU with multiple keys.

Positional arguments:
  n              degree of the underlying ring
  q              prime modulus

Optional arguments:
  -h, --help     shows help message and exits
  -v, --version  prints version information and exits
  -c, --coeffs   number of coefficients. 2 for binary, 3 for ternary [nargs=0..1] [default: 2]
  -s, --seed     optionally fix seed. If seed is -1 then use a random seed. [nargs=0..1] [default: -1]
  -r, --ring     use 1 for NTRU: x^n - 1, 2 for NTRU2: x^n + 1, 3 for NTRUPrime: x^n - x - 1 or 4 for NTTRU: x^n - x^(n/2) + 1. [nargs=0..1] [default: 1]
  --verbose      increase output verbosity

Subcommands:
  all           All-in-one command.
  keygen        Generate NTRU keys with a shared denominator.
  recover       Recover key from linearized system.
  system        Create linearized system.
  verify        Verify the files contain secret keys which are rotations of each other.
```

The `arora-ge-ntru` binary can only be run using one of the five subcommands
above.
For example, to generate 10 samples in the ring `Z[x]/(x^n - 1)` with `n = 16` and modulo `q = 31`
you can do:
```
./arora-ge-ntru 16 31 -c 2 -r 1 --verbose keygen -k 10 --pk_output=pk --sk_output=sk
```
This saves the public keys to `pk` and the secret denominator to `sk`.

The other subcommands are similar. For help with a subcommand, e.g. `recover`, do:
```
$ ./arora-ge-ntru 16 31 recover --help
Usage: recover [--help] [--version] --input VAR [--output VAR] [--nullonly]

Recover key from linearized system.

Optional arguments:
  -h, --help     shows help message and exits
  -v, --version  prints version information and exits
  -i, --input    input file of linearized system (output of system subcommand) [required]
  -o, --output   optional output file
  --nullonly     flag -- only output nullspace and then stop
```
(Note some argument for `n` and `q` is still required.)

# License
Licensed under the MIT License <http://opensource.org/licenses/MIT>.

This software uses the Argparse header-only library for argument parsing, licensed
under the MIT License and available at <http://github.com/p-ranav/argparse>.
