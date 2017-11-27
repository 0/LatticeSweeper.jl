module LatticeSweeper

using LinearMaps
using OffsetArrays
using TensorOperations

export
    MPS,
    MPO,
    SweepSchedule,
    dmrg!

"""
Sweep or contraction direction.
"""
@enum Direction Left Right

include("util.jl")

include("mps.jl")
include("mpo.jl")
include("contraction.jl")

include("settings.jl")

include("dmrg.jl")

end
