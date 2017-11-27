"""
    dmrg_step!{L,T}(
        H::MPO{L,T},
        psi::MPS{L,T},
        H_cntrctns::OffsetArray{Contraction{T,1},1},
        dir::Direction,
        site::Int)

Update `psi` by doing a two-site DMRG step in the direction `dir` at sites
`site` and `site+1` using the Hamiltonian `H` and Hamiltonian contractions
`H_cntrctns`.
"""
function dmrg_step!{L,T}(
        H::MPO{L,T},
        psi::MPS{L,T},
        H_cntrctns::OffsetArray{Contraction{T,1},1},
        dir::Direction,
        site::Int)

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
    # after a single iteration and do it twice.
    eigen = eigs(H_map; nev=1, which=:SR, tol=Inf, maxiter=1, v0=vec(prev))
    eigen = eigs(H_map; nev=1, which=:SR, tol=Inf, maxiter=1, v0=vec(eigen[2]))

    # Update the MPS.
    wf_mat = reshape(eigen[2], prod(size(prev, 1, 2)), prod(size(prev, 3, 4)))
    U, S, V = svd(wf_mat)
    if dir == Right
        A = U
        B = diagm(S) * V'
    elseif dir == Left
        A = U * diagm(S)
        B = V'
    end
    psi.tnsrs[site] = reshape(A, size(prev, 1), size(prev, 2), length(S))
    psi.tnsrs[site+1] = reshape(B, length(S), size(prev, 3), size(prev, 4))

    # Update the H contraction.
    if dir == Right
        H_cntrctns[site] = contract_site(H_cntrctns[site-1], psi, H)
    elseif dir == Left
        H_cntrctns[site+1] = contract_site(H_cntrctns[site+2], psi, H)
    end

    nothing
end


"""
    dmrg!{L,T}(psi::MPS{L,T}, H::MPO{L,T}, num_sweeps::Int)

Perform `num_sweeps` sweeps of two-site DMRG on the state `psi` using the
Hamiltonian `H`.
"""
function dmrg!{L,T}(psi::MPS{L,T}, H::MPO{L,T}, num_sweeps::Int)
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

    for n in 1:num_sweeps
        if n % 2 != 0
            # Right sweep.
            for i in 1:(L-2)
                dmrg_step!(H, psi, H_cntrctns, Right, i)
            end
        else
            # Left sweep.
            for i in (L-1):-1:2
                dmrg_step!(H, psi, H_cntrctns, Left, i)
            end
        end
    end

    if num_sweeps % 2 == 0
        # Ended on the left.
        H_cntrctn = contract_site(H_cntrctns[3], psi, H)
        E0 = contract_sides(H_cntrctns[1], H_cntrctn)
    else
        # Ended on the right.
        H_cntrctn = contract_site(H_cntrctns[L-2], psi, H)
        E0 = contract_sides(H_cntrctn, H_cntrctns[L])
    end

    realize(E0)
end
