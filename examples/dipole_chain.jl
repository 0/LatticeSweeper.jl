#!/usr/bin/env julia

# Dipole chain DMRG example.

push!(LOAD_PATH, joinpath(dirname(@__FILE__), "../src"))
using LatticeSweeper

using ArgParse

s = ArgParseSettings()
s.autofix_names = true
@add_arg_table s begin
    "-R"
        metavar = "R"
        help = "separation distance"
        arg_type = Float64
        required = true
    "-N"
        metavar = "N"
        help = "number of rotors"
        arg_type = Int
        required = true
    "--l-max"
        metavar = "L"
        help = "local basis truncation"
        arg_type = Int
        required = true
    "--max-sweeps"
        metavar = "S"
        help = "maximum number of sweeps"
        arg_type = Int
        required = true
end
c = parse_args(ARGS, s, as_symbols=true)

R = c[:R]
N = c[:N]
l_max = c[:l_max]
max_sweeps = c[:max_sweeps]

# l m = 0 0, 1 -1, 1 0, 1 1, ...
basis = [(l, m) for l in 0:l_max for m in -l:l]
basis_size = length(basis)

# One-site operators.
id = zeros(basis_size, basis_size)
L2 = zeros(basis_size, basis_size)
Lp = zeros(basis_size, basis_size)
Lm = zeros(basis_size, basis_size)
LpMp = zeros(basis_size, basis_size)
LpMm = zeros(basis_size, basis_size)
LmMp = zeros(basis_size, basis_size)
LmMm = zeros(basis_size, basis_size)
for (col, (l, m)) in enumerate(basis)
    for (row, (lp, mp)) in enumerate(basis)
        if lp == l && mp == m
            id[row, col] = 1.0
            L2[row, col] = l*(l+1)
        end
        if lp == l+1 && mp == m
            Lp[row, col] = sqrt((l-m+1)*(l+m+1)/((2l+1)*(2l+3)))
        end
        if lp == l-1 && mp == m
            Lm[row, col] = sqrt((l-m)*(l+m)/((2l-1)*(2l+1)))
        end
        if lp == l+1 && mp == m+1
            LpMp[row, col] = -sqrt((l+m+1)*(l+m+2)/((2l+1)*(2l+3)))
        end
        if lp == l+1 && mp == m-1
            LpMm[row, col] = -sqrt((l-m+1)*(l-m+2)/((2l+1)*(2l+3)))
        end
        if lp == l-1 && mp == m+1
            LmMp[row, col] = sqrt((l-m-1)*(l-m)/((2l-1)*(2l+1)))
        end
        if lp == l-1 && mp == m-1
            LmMm[row, col] = sqrt((l+m-1)*(l+m)/((2l-1)*(2l+1)))
        end
    end
end

# Dipole-dipole potential operator.
ks = Float64[]
op1s = Matrix{Float64}[]
op2s = Matrix{Float64}[]
for (k, OP1, OP2) in [(-2.0, [Lp, Lm], [Lp, Lm]),
                      (-0.5, [LpMm, LmMm], [LpMp, LmMp]),
                      (-0.5, [LpMp, LmMp], [LpMm, LmMm])]
    for op1 in OP1
        for op2 in OP2
            push!(ks, k)
            push!(op1s, op1)
            push!(op2s, op2)
        end
    end
end
K = length(ks)

# Dipole chain Hamiltonian.
H_tnsr = zeros(2+(N-1)*K, 2+(N-1)*K, basis_size, basis_size)
H_tnsr[1, 1, :, :] = id
H_tnsr[2, 1, :, :] = L2
H_tnsr[2, 2, :, :] = id
for i in 1:(N-1)
    for k in 1:K
        H_tnsr[2, 2+(i-1)*K+k, :, :] = ks[k]*op1s[k]/(R*i)^3
    end
end
for k in 1:K
    H_tnsr[2+k, 1, :, :] = op2s[k]
end
for i in 1:(N-2)
    for k in 1:K
        H_tnsr[2+i*K+k, 2+(i-1)*K+k, :, :] = id
    end
end
H = MPO(H_tnsr, N)
compress!(H, 1e-15)

# Starting wavefunction.
wf = zeros(basis_size)
wf[1] = 1.0
psi = MPS(wf, N)

hist = dmrg!(psi, H, SweepSchedule(max_sweeps))
hist.converged || warn("Ground state not converged.")

# Ground state energy.
println("E0 = $(hist[end].energy)")
# Von Neumann entanglement entropy.
println("SvN = $(S_vn(hist[end].middle_eigvals))")
