"""
    SweepDetails

Details of a single DMRG sweep.
"""
struct SweepDetails
    "Total energy of the state."
    energy::Float64
    "Normalized fluctuation of the energy: `<(H - <H>)^2>/<H>^2`."
    dH2::Float64
    "Eigenvalues across the middle bond."
    middle_eigvals::Vector{Float64}

    "Duration of the sweep in seconds."
    duration::Float64
end

"""
    SweepDetails(state::SweepState)

Create a `SweepDetails` using the values in `state`.
"""
SweepDetails(state::SweepState) =
    SweepDetails(state.energy, state.dH2, state.middle_eigvals, state.duration)


"""
    SweepHistory

Information about a sequence DMRG sweeps.
"""
struct SweepHistory
    "Details from all completed sweeps."
    details::Vector{SweepDetails}
    "Whether the state is converged to the desired tolerance."
    converged::Bool
end

Base.endof(hist::SweepHistory) = length(hist.details)
Base.getindex(hist::SweepHistory, i::Int) = hist.details[i]
Base.length(hist::SweepHistory) = length(hist.details)
