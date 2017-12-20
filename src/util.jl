"""
    realize(x::Number)

Make real, but carefully.
"""
realize(x::Number) = error("Unrecognized numerical type")
realize(x::Real) = x
function realize(x::Complex)
    abs(imag(x)) < 1e-14 || error("Non-zero imaginary component")
    real(x)
end


"""
    S_vn(eigvals::Vector{Float64})

Von Neumann entropy of `eigvals`.
"""
S_vn(eigvals::Vector{Float64}) = -sum(x * log(x) for x in eigvals if x > 0)
