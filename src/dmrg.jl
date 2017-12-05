"""
    dmrg_step!{L,T}(state::SweepState{L,T}, set::SweepSettings)

Update `state` by doing a two-site DMRG step with the parameters in `set`.
"""
function dmrg_step!{L,T}(state::SweepState{L,T}, set::SweepSettings)
    psi, H, site = state.psi, state.H, state.site

    left = state.H_cntrctns[site-1].tnsr
    right = state.H_cntrctns[site+2].tnsr

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
    # Reduced density matrix eigenvalues.
    eigvals = S[1:trunc_len].^2

    # Update the MPS.
    if state.dir == Right
        A = U[:, 1:trunc_len]
        B = diagm(S[1:trunc_len]) * V[:, 1:trunc_len]'
    elseif state.dir == Left
        A = U[:, 1:trunc_len] * diagm(S[1:trunc_len])
        B = V[:, 1:trunc_len]'
    end
    psi.tnsrs[site] = reshape(A, size(prev, 1, 2)..., trunc_len)
    psi.tnsrs[site+1] = reshape(B, trunc_len, size(prev, 3, 4)...)

    # Update the contractions.
    if state.dir == Right
        state.H_cntrctns[site] = state.H_cntrctns[site-1] * (psi, H)
    elseif state.dir == Left
        state.H_cntrctns[site+1] = state.H_cntrctns[site+2] * (psi, H)
    end
    state.H2_cntrctn = state.H2_cntrctn * (psi, H, H)

    # Update the maximum bond dimension.
    state.mps_max_bond_dim = maximum(size(tnsr, 3) for tnsr in psi.tnsrs)

    eigvals
end


"""
    dmrg!{L,T}(psi::MPS{L,T}, H::MPO{L,T}, sch::SweepSchedule;
               outputs::Vector{<:SweepOutput}=SweepOutput[])

Perform two-site DMRG on the state `psi` using the Hamiltonian `H` with
parameters specified in `sch`.
"""
function dmrg!{L,T}(psi::MPS{L,T}, H::MPO{L,T}, sch::SweepSchedule;
                    outputs::Vector{<:SweepOutput}=SweepOutput[])
    # At least 2 sites.
    L >= 2 || throw(DomainError())

    state = SweepState(psi, H)
    sweep_details = SweepDetails[]

    init.(outputs)

    converged = false
    for n in 1:sch.max_sweeps
        sweep_start_time = time()

        set = sch[n]

        state.sweep_num = n
        if L == 2
            state.dir = Right
        else
            state.dir = flip(state.dir)
        end
        state.H2_cntrctn = cap_contraction(T, state.dir, L, 2)
        if L == 2
            # 1 2
            # ^ ^
            range = 1:1
        else
            if state.dir == Right
                # 1 2 3 4 5
                # ^ ^
                #   ^ ^
                #     ^ ^
                range = 1:(L-2)
            elseif state.dir == Left
                # 1 2 3 4 5
                #       ^ ^
                #     ^ ^
                #   ^ ^
                range = (L-1):-1:2
            end
        end

        # Sweep!
        for i in range
            state.site = i
            eigvals = dmrg_step!(state, set)
            i == div(L, 2) && (state.middle_eigvals = eigvals)

            step.(outputs, state)
        end

        if L == 2
            energy = state.H_cntrctns[1] * (psi, H) * state.H_cntrctns[L+1]
            H2 = state.H2_cntrctn * (psi, H, H) *
                 cap_contraction(T, Left, L, 2)
        else
            if state.dir == Right
                energy = state.H_cntrctns[L-2] * (psi, H) * state.H_cntrctns[L]
                H2 = state.H2_cntrctn * (psi, H, H) * (psi, H, H) *
                     cap_contraction(T, Left, L, 2)
            elseif state.dir == Left
                energy = state.H_cntrctns[1] * (psi, H) * state.H_cntrctns[3]
                H2 = cap_contraction(T, Right, L, 2) *
                     (psi, H, H) * (psi, H, H) * state.H2_cntrctn
            end
        end

        state.energy = realize(energy)
        state.dH2 = realize(H2)/energy^2 - 1.0
        state.duration = time() - sweep_start_time

        push!(sweep_details, SweepDetails(set, state))

        sweep.(outputs, state)

        if abs(state.dH2) <= sch.tolerance
            converged = true
            break
        end
    end

    done.(outputs)

    SweepHistory(sweep_details, converged)
end
