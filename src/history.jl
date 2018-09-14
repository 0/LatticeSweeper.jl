"""
    SweepDetails

Details of a single DMRG sweep.
"""
struct SweepDetails
    "Settings for the sweep."
    set::SweepSettings

    "Total energy of the state."
    energy::Float64
    "Normalized fluctuation of the energy: `<(H - <H>)^2>/<H>^2`."
    dH2::Float64
    "Eigenvalues across the middle bond."
    middle_eigvals::Vector{Float64}

    "Maximum wavefunction bond dimension."
    mps_max_bond_dim::Int

    "Duration of the sweep in seconds."
    duration::Float64
end

"""
    SweepDetails(set::SweepSettings, state::SweepState)

Create a `SweepDetails` for settings `set` using the values in `state`.
"""
SweepDetails(set::SweepSettings, state::SweepState) =
    SweepDetails(set, state.energy, state.dH2, state.middle_eigvals,
                 state.mps_max_bond_dim, state.duration)


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

Base.getindex(hist::SweepHistory, i::Int) = hist.details[i]
Base.lastindex(hist::SweepHistory) = length(hist.details)
Base.length(hist::SweepHistory) = length(hist.details)
