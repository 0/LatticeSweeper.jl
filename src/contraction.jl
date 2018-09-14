"""
    Contraction{T<:Number,N,dir}

Contraction of 2 `MPS`s and `N` `MPO`s from an edge, going in the direction
`dir` with elements of type `T`.

For example, a `Contraction` to the right with 1 MPO over 2 sites would look
like

          cap   1    2
    #-     #-  -#-  -#-
    #      #    |    |
    #      #
    #      #    |    |
    #-  =  #-  -#-  -#-
    #      #    |    |
    #      #
    #      #    |    |
    #-     #-  -#-  -#-
"""
struct Contraction{T<:Number,N,dir}
    "Last site added to the contraction."
    site::Int
    "Contracted tensor."
    tnsr::Array{T}
end

"""
    Contraction(dir::Direction, site::Int, tnsr::Array)

Make a `Contraction` up to site `site` going in the direction `dir` and
containing the elements `tnsr`.
"""
Contraction(dir::Direction, site::Int, tnsr::Array{T,N}) where {T,N} =
    Contraction{T,N-2,dir}(site, tnsr)

"""
    cap_contraction(::Type{T}, dir::Direction, L::Int, num_mpos::Int)

Make a cap `Contraction` of type `T` going in the direction `dir` for `L` sites
and `num_mpos` `MPO`s.
"""
cap_contraction(::Type{T}, dir::Direction, L::Int, num_mpos::Int) where {T} =
    Contraction(dir, dir == Right ? 0 : L+1, ones(T, ones(Int, num_mpos+2)...))

"""
    dummy_contraction(::Type{T}, num_mpos::Int)

Make a dummy `Contraction` of type `T` for `num_mpos` `MPO`s.
"""
dummy_contraction(::Type{T}, num_mpos::Int) where {T} =
    Contraction(Left, -1, zeros(T, ones(Int, num_mpos+2)...))


"""
    contract_site(
        cntrctn::Contraction{T,N,dir},
        s1::MPS{L,T},
        s2::MPS{L,T},
        os::Vararg{MPO{L,T},N})

Contract the next site of `MPS` `s1`, `MPS` `s2`, and `MPO`s `os` into
`cntrctn`, producing a new `Contraction`.
"""
@generated function contract_site(
        cntrctn::Contraction{T,N,dir},
        s1::MPS{L,T},
        s2::MPS{L,T},
        os::Vararg{MPO{L,T},N}) where {T,N,dir,L}

    # tnsr[i_s1, i_1, ..., i_N, i_s2]
    input = :(tnsr[i_s1])
    # result[j_s1, j_1, ..., j_N, j_s2]
    output = :(result[j_s1])
    for idx = 1:N
        push!(input.args, Symbol("i_", idx))
        push!(output.args, Symbol("j_", idx))
    end
    push!(input.args, :i_s2)
    push!(output.args, :j_s2)

    label_left, label_right = dir == Right ? ("i_", "j_") : ("j_", "i_")

    # All the tensors multiplied on the right hand side.
    rhs = Expr(:call, :*, input)
    left = Symbol(label_left, "s1")
    nxt = Symbol("k_", 0)
    right = Symbol(label_right, "s1")
    push!(rhs.args, :(s1.tnsrs[site][$left, $nxt, $right]))
    for idx = 1:N
        left = Symbol(label_left, idx)
        prv = Symbol("k_", idx-1)
        nxt = Symbol("k_", idx)
        right = Symbol(label_right, idx)
        push!(rhs.args, :(os[$idx].tnsrs[site][$left, $prv, $nxt, $right]))
    end
    left = Symbol(label_left, "s2")
    prv = Symbol("k_", N)
    right = Symbol(label_right, "s2")
    if T <: Real
        push!(rhs.args, :(s2.tnsrs[site][$left, $prv, $right]))
    elseif T <: Complex
        push!(rhs.args, :(conj(s2.tnsrs[site][$left, $prv, $right])))
    end

    quote
        tnsr = cntrctn.tnsr
        # Next site index.
        site = cntrctn.site + (dir == Right ? 1 : -1)

        # When dir == Right and T == Float64:
        #   result[j_s1, j_1, ..., j_N, j_s2]
        #     := tnsr[i_s1, i_1, ..., i_N, i_s2]
        #        * s1.tnsrs[site][i_s1, k_0, j_s1]
        #        * os[1].tnsrs[site][i_1, k_0, k_1, j_1]
        #        * ...
        #        * os[N].tnsrs[site][i_N, k_(N-1), k_N, j_N]
        #        * s2.tnsrs[site][i_s2, k_N, j_s2]
        @tensor $output := $rhs

        Contraction(dir, site, result)
    end
end

"""
    contract_site(
        cntrctn::Contraction{T,N,dir},
        s::MPS{L,T},
        os::Vararg{MPO{L,T},N})

Contract the next site of `MPS` `s` (on both sides) and `MPO`s `os` into
`cntrctn`, producing a new `Contraction`.
"""
contract_site(
        cntrctn::Contraction{T,N,dir},
        s::MPS{L,T},
        os::Vararg{MPO{L,T},N}) where {T,N,dir,L} =
    contract_site(cntrctn, s, s, os...)

Base.:*(cntrctn::Contraction, args::Tuple{MPS,Vararg{MPO}}) =
    contract_site(cntrctn, args...)


"""
    contract_sides(
        cntrctn_left::Contraction{T,N,Right},
        cntrctn_right::Contraction{T,N,Left})

Contract together left and right sides, producing a scalar.
"""
@generated function contract_sides(
        cntrctn_left::Contraction{T,N,Right},
        cntrctn_right::Contraction{T,N,Left}) where {T,N}

    # tnsr_left[i_s1, i_1, ..., i_N, i_s2]
    left = :(tnsr_left[i_s1])
    # tnsr_right[i_s1, i_1, ..., i_N, i_s2]
    right = :(tnsr_right[i_s1])
    for idx = 1:N
        push!(left.args, Symbol("i_", idx))
        push!(right.args, Symbol("i_", idx))
    end
    push!(left.args, :i_s2)
    push!(right.args, :i_s2)

    quote
        @assert cntrctn_left.site == cntrctn_right.site - 1

        tnsr_left = cntrctn_left.tnsr
        tnsr_right = cntrctn_right.tnsr

        # tnsr_left[i_s1, i_1, ..., i_N, i_s2]
        # * tnsr_right[i_s1, i_1, ..., i_N, i_s2]
        @tensor $left * $right
    end
end

Base.:*(cntrctn_left::Contraction, cntrctn_right::Contraction) =
    contract_sides(cntrctn_left, cntrctn_right)
