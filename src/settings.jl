"""
    ForeverVector{T}

Vector of values of type `T` that reuses the value at the last index for
subsequent indices.
"""
struct ForeverVector{T}
    "Elements of the vector."
    xs::Vector{T}
end

"""
    ForeverVector(x)

Create a `ForeverVector` with value(s) `x`.
"""
ForeverVector(x) = ForeverVector([x])

Base.getindex(fv::ForeverVector, i::Int) =
    i <= length(fv.xs) ? fv.xs[i] : fv.xs[end]


"""
    SweepSettings

Parameters for a single DMRG sweep.
"""
struct SweepSettings
    "Number of diagonalization iterations in a sweep."
    num_iters::Int
end


"""
    SweepSchedule

Schedule of parameters for multiple DMRG sweeps.
"""
struct SweepSchedule
    "Number of sweeps."
    num_sweeps::Int

    # SweepSettings parameters.
    num_iters::ForeverVector{Int}
end

"""
    SweepSchedule(num_sweeps; kwargs...)

Create a `SweepSchedule` for `num_sweeps` sweeps, with optional parameters
specified in `kwargs`.
"""
function SweepSchedule(num_sweeps; kwargs...)
    # At least 1 sweep.
    num_sweeps >= 1 || throw(DomainError())

    kwargs = Dict(kwargs)

    num_iters = ForeverVector(get(kwargs, :num_iters, 2))

    SweepSchedule(num_sweeps, num_iters)
end

Base.getindex(sch::SweepSchedule, i::Int) = SweepSettings(sch.num_iters[i])
