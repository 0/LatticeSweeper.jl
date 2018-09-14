"""
    MPS{L,T<:Number}

Matrix product state for `L` sites with numerical entries of type `T`.

The tensors are stored with the indices `[a, m, b]` corresponding to

       m
       |
    a -#- b

where the physical index (`m`) is between the bond indices (`a`, `b`). The end
tensors also have this structure, so they must be capped explicitly.
"""
struct MPS{L,T<:Number}
    "Tensors making up the MPS."
    tnsrs::Vector{Array{T,3}}
end

"""
    MPS(vs::Vector{Vector{T}})

Create an `MPS` representing a product state (all bonds have dimension 1),
where each site is described by the corresponding element of `vs`.
"""
function MPS(vs::Vector{Vector{T}}) where {T}
    L = length(vs)

    tnsrs = Vector{Array{T,3}}(undef, L)
    for i in 1:L
        tnsrs[i] = reshape(copy(vs[i]), 1, :, 1)
    end

    MPS{L,T}(tnsrs)
end

"""
    MPS(v::Vector{T}, L)

Create an `MPS` for `L` sites representing a uniform product state (all bonds
have dimension 1), where each site is described by `v`.
"""
MPS(v::Vector{T}, L) where {T} = MPS([v for _ in 1:L])
