"""
    dmrg_step!{L,T}(
        H::MPO{L,T},
        psi::MPS{L,T},
        H_cntrctns::OffsetArray{Contraction{T,1},1},
        dir::Direction,
        site::Int,
        set::SweepSettings)

Update `psi` by doing a two-site DMRG step in the direction `dir` at sites
`site` and `site+1` using the Hamiltonian `H` and Hamiltonian contractions
`H_cntrctns`.
"""
function dmrg_step!{L,T}(
        H::MPO{L,T},
        psi::MPS{L,T},
        H_cntrctns::OffsetArray{Contraction{T,1},1},
        dir::Direction,
        site::Int,
        set::SweepSettings)

    left = H_cntrctns[site-1].tnsr
    right = H_cntrctns[site+2].tnsr

    # The previous two-site wavefunction will be used as the starting vector:
    #
    #        m     n
    #        |     |
    #     a -#- b -#- c
    @tensor prev[a, m, n, c] :=
        psi.tnsrs[site][a, m, b] * psi.tnsrs[site+1][b, n, c]

    # Application of the two-site-projected Hamiltonian onto a vector:
    #
    #     #- c           i -#
    #     #     n     p     #
    #     #     |     |     #
    #     #- b -#- d -#- h -#
    #     #     |     |     #
    #     #     m     o     #
    #     #     |     |     #
    #     #- a -#######- g -#
    function H_f(dest::AbstractVector, src::AbstractVector)
        dest_ = reshape(dest, size(prev)...)
        src_ = reshape(src, size(prev)...)
        @tensor dest_[c, n, p, i] =
            left[a, b, c] * src_[a, m, o, g] *
            H.tnsrs[site][b, m, n, d] * H.tnsrs[site+1][d, o, p, h] *
            right[g, h, i]
        dest
    end
    H_map = LinearMap{T}(H_f, prod(size(prev)); issymmetric=true)

    # Improve our guess for the ground state of the two-site-projected
    # Hamiltonian by using a sparse diagonalizer. Since there's no way to
    # request a specific number of iterations, we just say we're converged
    # after a single iteration and repeat several times.
    eigen = eigs(H_map; nev=1, which=:SR, tol=Inf, maxiter=1, v0=vec(prev))
    wf = vec(eigen[2])
    for _ in 1:(set.num_iters-1)
        eigen = eigs(H_map; nev=1, which=:SR, tol=Inf, maxiter=1, v0=wf)
        wf = vec(eigen[2])
    end

    # Update the MPS.
    wf_mat = reshape(wf, prod(size(prev, 1, 2)), prod(size(prev, 3, 4)))
    U, S, V = svd(wf_mat)
    trunc_len = length(S)
    trunc_cutoff = 0.0
    while trunc_len > set.bond_max ||
          (set.bond_min < trunc_len &&
           trunc_cutoff + S[trunc_len]^2 <= set.cutoff_max)

        trunc_cutoff += S[trunc_len]^2
        trunc_len -= 1
    end
    # Because the diagonalization results in a normalized eigenvector, we are
    # guaranteed that sum(S.^2) == 1. However, when we throw away some singular
    # values, we should renormalize.
    S ./= sqrt(sum(S[1:trunc_len].^2))
    if dir == Right
        A = U[:, 1:trunc_len]
        B = diagm(S[1:trunc_len]) * V[:, 1:trunc_len]'
    elseif dir == Left
        A = U[:, 1:trunc_len] * diagm(S[1:trunc_len])
        B = V[:, 1:trunc_len]'
    end
    psi.tnsrs[site] = reshape(A, size(prev, 1), size(prev, 2), trunc_len)
    psi.tnsrs[site+1] = reshape(B, trunc_len, size(prev, 3), size(prev, 4))

    # Update the H contraction.
    if dir == Right
        H_cntrctns[site] = contract_site(H_cntrctns[site-1], psi, H)
    elseif dir == Left
        H_cntrctns[site+1] = contract_site(H_cntrctns[site+2], psi, H)
    end

    # Reduced density matrix eigenvalues.
    S[1:trunc_len].^2
end


"""
    dmrg!{L,T}(psi::MPS{L,T}, H::MPO{L,T}, sch::SweepSchedule)

Perform two-site DMRG on the state `psi` using the Hamiltonian `H` with
parameters specified in `sch`.
"""
function dmrg!{L,T}(psi::MPS{L,T}, H::MPO{L,T}, sch::SweepSchedule)
    # At least 3 sites.
    L >= 3 || throw(DomainError())

    # Initialize the H contractions for the first sweep.
    H_cntrctns = OffsetArray(Contraction{T,1}, 0:(L+1))
    H_cntrctns[L+1] = cap_contraction(T, Left, L, 1)
    for i in L:-1:3
        H_cntrctns[i] = contract_site(H_cntrctns[i+1], psi, H)
    end
    H_cntrctns[2] = dummy_contraction(T, 1)
    H_cntrctns[1] = dummy_contraction(T, 1)
    H_cntrctns[0] = cap_contraction(T, Right, L, 1)

    sweep_details = SweepDetails[]
    middle_eigvals = Float64[]

    for n in 1:sch.num_sweeps
        if n % 2 != 0
            # Right sweep.
            for i in 1:(L-2)
                eigvals = dmrg_step!(H, psi, H_cntrctns, Right, i, sch[n])
                i == div(L, 2) && (middle_eigvals = eigvals)
            end
            energy = contract_site(H_cntrctns[L-2], psi, H) * H_cntrctns[L]
        else
            # Left sweep.
            for i in (L-1):-1:2
                eigvals = dmrg_step!(H, psi, H_cntrctns, Left, i, sch[n])
                i == div(L, 2) && (middle_eigvals = eigvals)
            end
            energy = H_cntrctns[1] * contract_site(H_cntrctns[3], psi, H)
        end

        push!(sweep_details, SweepDetails(realize(energy), middle_eigvals))
    end

    SweepHistory(sweep_details)
end
