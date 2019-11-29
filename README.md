# LatticeSweeper

Simple two-site DMRG.
No support for quantum number conservation.

Tested with Julia 1.3.


## Installation

```
pkg> add https://github.com/0/LatticeSweeper.jl.git
```

In order to run the example scripts in `examples/`, you will also need to
```
pkg> add ArgParse
```

### Application project

If you're working with a clone of this repository, you can use the basic application project in `examples/`, which already has both `LatticeSweeper` and `ArgParse` as dependencies.
From the repository root, run
```
julia --project=examples
```
and then
```
pkg> dev .
```
to create `examples/Manifest.toml` with a development version of `LatticeSweeper`.


## Examples

To run the following examples, you should set the project (e.g. using `--project` or `JULIA_PROJECT`) to a Julia project that has the prerequisites installed.

* `julia examples/tfim.jl --help`
* `julia --color=yes examples/tfim.jl -g 1.0 -L 8 --max-sweeps 4`
* `julia examples/dipole_chain.jl --help`
* `julia --color=yes examples/dipole_chain.jl -R 1.0 -N 3 --l-max 1 --max-sweeps 4`


## References

The DMRG implementation is based largely on [Schollwöck's 2011 review "The density-matrix renormalization group in the age of matrix product states"](http://www.sciencedirect.com/science/article/pii/S0003491610001752) ([preprint](https://arxiv.org/abs/1008.3477v2)).


## License

Provided under the terms of the MIT license.
See `LICENSE` for more information.
