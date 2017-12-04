# LatticeSweeper.jl

Two-site DMRG without any bells or whistles.

Tested with Julia 0.6.


## Installation

1. `Pkg.clone("https://github.com/0/LatticeSweeper.jl.git")`


### Requirements

These should be pulled in automatically when installing this package.
To use it without installing it (e.g. from a local git checkout), you'll need to manually obtain the following dependencies:

* ArgParse (`Pkg.add("ArgParse")`)
* LinearMaps (`Pkg.add("LinearMaps")`)
* OffsetArrays (`Pkg.add("OffsetArrays")`)
* TensorOperations (`Pkg.add("TensorOperations")`)


## Examples

* `julia examples/tfim.jl --help`
* `julia --color=yes examples/tfim.jl -g 1.0 -L 8 --max-sweeps 4`
* `julia examples/dipole_chain.jl --help`
* `julia --color=yes examples/dipole_chain.jl -R 1.0 -N 3 --l-max 1 --max-sweeps 4`


## References

The DMRG implementation is based largely on [Schollwöck's 2011 review "The density-matrix renormalization group in the age of matrix product states"](http://www.sciencedirect.com/science/article/pii/S0003491610001752) ([preprint](https://arxiv.org/abs/1008.3477v2)).


## License

Provided under the terms of the MIT license.
See `LICENSE` for more information.
