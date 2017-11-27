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
