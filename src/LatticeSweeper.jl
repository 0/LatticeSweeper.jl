module LatticeSweeper

using LinearMaps
using OffsetArrays
using TensorOperations

export
    MPS,
    MPO,
    dmrg!

"""
Sweep or contraction direction.
"""
@enum Direction Left Right

include("util.jl")

include("mps.jl")
include("mpo.jl")
include("contraction.jl")

include("dmrg.jl")

end
