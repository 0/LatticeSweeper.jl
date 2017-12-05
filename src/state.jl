"""
    SweepState

State of a DMRG sweep.
"""
mutable struct SweepState{L,T}
    "Wavefunction."
    psi::MPS{L,T}
    "Hamiltonian."
    H::MPO{L,T}

    "Site index."
    site::Int
    "Sweep direction."
    dir::Direction

    "Hamiltonian contractions."
    H_cntrctns::OffsetArray{Contraction{T,1},1}
    "Hamiltonian fluctuation contraction."
    H2_cntrctn::Contraction{T,2}

    # SweepDetails parameters.
    energy::Float64
    dH2::Float64
    middle_eigvals::Vector{Float64}
    duration::Float64
end

"""
    SweepState{L,T}(psi::MPS{L,T}, H::MPO{L,T})

Generate a blank `SweepState` for wavefunction `psi` and Hamiltonian `H`.
"""
function SweepState{L,T}(psi::MPS{L,T}, H::MPO{L,T})
    # Initialize the H contractions for the first sweep.
    H_cntrctns = OffsetArray(Contraction{T,1}, 0:(L+1))
    H_cntrctns[L+1] = cap_contraction(T, Left, L, 1)
    for i in L:-1:3
        H_cntrctns[i] = H_cntrctns[i+1] * (psi, H)
    end
    H_cntrctns[2] = dummy_contraction(T, 1)
    H_cntrctns[1] = dummy_contraction(T, 1)
    H_cntrctns[0] = cap_contraction(T, Right, L, 1)

    # Use dummy values for parameters that we don't know yet.
    SweepState(psi, H, -1, Left, H_cntrctns, dummy_contraction(T, 2),
               NaN, NaN, Float64[], NaN)
end
