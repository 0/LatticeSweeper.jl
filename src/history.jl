"""
    SweepDetails

Details of a single DMRG sweep.
"""
struct SweepDetails
    "Total energy of the state."
    energy::Float64
    "Eigenvalues across the middle bond."
    middle_eigvals::Vector{Float64}
end


"""
    SweepHistory

Information about a sequence DMRG sweeps.
"""
struct SweepHistory
    "Details from all completed sweeps."
    details::Vector{SweepDetails}
end

Base.endof(hist::SweepHistory) = length(hist.details)
Base.getindex(hist::SweepHistory, i::Int) = hist.details[i]
Base.length(hist::SweepHistory) = length(hist.details)
