#!/usr/bin/env julia

# Transverse field Ising model DMRG example.

push!(LOAD_PATH, joinpath(dirname(@__FILE__), "../src"))
using LatticeSweeper

using ArgParse

s = ArgParseSettings()
s.autofix_names = true
@add_arg_table s begin
    "-g"
        metavar = "g"
        help = "transverse field strength"
        arg_type = Float64
        required = true
    "-L"
        metavar = "L"
        help = "number of sites"
        arg_type = Int
        required = true
    "--num-sweeps"
        metavar = "S"
        help = "number of sweeps"
        arg_type = Int
        required = true
end
c = parse_args(ARGS, s, as_symbols=true)

g = c[:g]
L = c[:L]
num_sweeps = c[:num_sweeps]

# Pauli matrices.
id = [1.0 0.0; 0.0 1.0]
sigma_x = [0.0 1.0; 1.0 0.0]
sigma_z = [1.0 0.0; 0.0 -1.0]

# Transverse field Ising model Hamiltonian.
H_tnsr = zeros(3, 3, 2, 2)
H_tnsr[1, 1, :, :] = id
H_tnsr[2, 1, :, :] = -g*sigma_x
H_tnsr[2, 2, :, :] = id
H_tnsr[2, 3, :, :] = -sigma_z
H_tnsr[3, 1, :, :] = sigma_z
H = MPO(H_tnsr, L)

# Starting wavefunction.
spin_up = [1.0, 0.0]
psi = MPS(spin_up, L)

hist = dmrg!(psi, H, SweepSchedule(num_sweeps))

# Ground state energy.
println("E0 = $(hist[end].energy)")
# Von Neumann entanglement entropy.
println("SvN = $(S_vn(hist[end].middle_eigvals))")
