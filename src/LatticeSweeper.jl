module LatticeSweeper

using LinearMaps
using OffsetArrays
using TensorOperations

export
    S_vn,

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
include("history.jl")

include("dmrg.jl")

end
