"""
    Direction

Sweep or contraction direction.
"""
@enum Direction Left Right


"""
    flip(dir::Direction)

Flip `dir` (`Left` to `Right` and back).
"""
flip(dir::Direction) = dir == Left ? Right : Left


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


"""
    println_result(msg...)

Output a result so that it stands out.
"""
println_result(msg...) = print_with_color(:green, msg..., "\n", bold=true)
