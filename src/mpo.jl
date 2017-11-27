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
