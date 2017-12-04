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

    "Minimum bond dimension to keep when truncating."
    bond_min::Int
    "Maximum bond dimension to keep when truncating."
    bond_max::Int
    "Maximum eigenvalue weight to drop when truncating."
    cutoff_max::Float64
end


"""
    SweepSchedule

Schedule of parameters for multiple DMRG sweeps.
"""
struct SweepSchedule
    "Maximum number of sweeps."
    max_sweeps::Int
    "Wavefunction tolerance."
    tolerance::Float64

    # SweepSettings parameters.
    num_iters::ForeverVector{Int}
    bond_min::ForeverVector{Int}
    bond_max::ForeverVector{Int}
    cutoff_max::ForeverVector{Float64}
end

"""
    SweepSchedule(max_sweeps::Int; tolerance=nothing, kwargs...)

Create a `SweepSchedule` for `max_sweeps` sweeps, with optional parameters
specified in `kwargs`.
"""
function SweepSchedule(max_sweeps::Int; tolerance=nothing, kwargs...)
    # At least 1 sweep.
    max_sweeps >= 1 || throw(DomainError())

    if tolerance !== nothing
        # Non-negative tolerance.
        tolerance >= 0.0 || throw(DomainError())
    else
        # Default.
        tolerance = 1e-5
    end

    kwargs = Dict(kwargs)

    num_iters = ForeverVector(get(kwargs, :num_iters, 2))
    bond_min = ForeverVector(get(kwargs, :bond_min, 5))
    bond_max = ForeverVector(get(kwargs, :bond_max, 100))
    cutoff_max = ForeverVector(get(kwargs, :cutoff_max, 1e-10))

    SweepSchedule(max_sweeps, tolerance,
                  num_iters, bond_min, bond_max, cutoff_max)
end

Base.getindex(sch::SweepSchedule, i::Int) =
    SweepSettings(sch.num_iters[i],
                  sch.bond_min[i],
                  sch.bond_max[i],
                  sch.cutoff_max[i])
