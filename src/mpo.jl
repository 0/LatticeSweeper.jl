"""
    MPO{L,T<:Number}

Matrix product operator for `L` sites with numerical entries of type `T`.

The tensors are stored with the indices `[a, m, n, b]` corresponding to

       n
       |
    a -#- b
       |
       m

where the physical indices (`m`, `n`) are between the bond indices (`a`, `b`).
The end tensors also have this structure, so they must be capped explicitly.
"""
struct MPO{L,T<:Number}
    "Tensors making up the MPO."
    tnsrs::Vector{Array{T,4}}
end

"""
    MPO{T}(x::Array{T,4}, L)

Create an `MPO` for `L` sites with all interior sites containing the tensor
`x`. The tensor is assumed to have the usual matrix-of-operators structure,
with the first two indices being the bond (matrix) dimension and the last two
indices being the physical (operator) dimension. The first and last sites only
use the second (not last!) row and first column of `x`, respectively.
"""
function MPO{T}(x::Array{T,4}, L)
    # At least 2 sites.
    L >= 2 || throw(DomainError())

    tnsrs = Vector{Array{T,4}}(L)
    # Row vector.
    tnsrs[1] = permutedims(x[2:2, :, :, :], (1, 3, 4, 2))
    for i in 2:(L-1)
        # Matrix.
        tnsrs[i] = permutedims(x, (1, 3, 4, 2))
    end
    # Column vector.
    tnsrs[L] = permutedims(x[:, 1:1, :, :], (1, 3, 4, 2))

    MPO{L,T}(tnsrs)
end


"""
    compress!{L}(o::MPO{L}, cutoff_max::Float64)

Compress `o` by throwing away singular values below `cutoff_max` at each bond.
"""
function compress!{L}(o::MPO{L}, cutoff_max::Float64)
    for range in [1:(L-2), (L-1):-1:1]
        for site in range
            # Combine adjacent tensors:
            #
            #        n     p
            #        |     |
            #     a -#- b -#- c
            #        |     |
            #        m     o
            @tensor M[a, m, n, o, p, c] :=
                o.tnsrs[site][a, m, n, b] * o.tnsrs[site+1][b, o, p, c]
            M_mat = reshape(M, prod(size(M, 1, 2, 3)), prod(size(M, 4, 5, 6)))

            U, S, V = svd(M_mat)
            # Operators are not normalized in general, so we should not expect
            # that sum(S.^2) == 1.
            trunc_total = sum(S.^2)
            trunc_len = length(S)
            trunc_cutoff = 0.0
            while (1 < trunc_len &&
                   trunc_cutoff + S[trunc_len]^2/trunc_total <= cutoff_max)
                trunc_cutoff += S[trunc_len]^2/trunc_total
                trunc_len -= 1
            end
            A = U[:, 1:trunc_len]
            B = diagm(S[1:trunc_len]) * V[:, 1:trunc_len]'
            o.tnsrs[site] = reshape(A, size(M, 1, 2, 3)..., trunc_len)
            o.tnsrs[site+1] = reshape(B, trunc_len, size(M, 4, 5, 6)...)
        end
    end

    nothing
end
