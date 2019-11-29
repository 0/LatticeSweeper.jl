#!/usr/bin/env julia

# Transverse field Ising model DMRG example.

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
    "--max-sweeps"
        metavar = "S"
        help = "maximum number of sweeps"
        arg_type = Int
        required = true
    "--output-file"
        metavar = "PATH"
        help = "path to output file"
        arg_type = String
end
c = parse_args(ARGS, s, as_symbols=true)

g = c[:g]
L = c[:L]
max_sweeps = c[:max_sweeps]
output_file = c[:output_file]

# Pauli matrices.
id = [1.0 0.0; 0.0 1.0]
sigma_x = [0.0 1.0; 1.0 0.0]
sigma_z = [1.0 0.0; 0.0 -1.0]

println("[ ] Constructing Hamiltonian.")
H_tnsr = zeros(3, 3, 2, 2)
H_tnsr[1, 1, :, :] = id
H_tnsr[2, 1, :, :] = -g*sigma_x
H_tnsr[2, 2, :, :] = id
H_tnsr[2, 3, :, :] = -sigma_z
H_tnsr[3, 1, :, :] = sigma_z
@time H = MPO(H_tnsr, L)
println("[+] Constructed Hamiltonian.")

println("[ ] Constructing wavefunction.")
spin_up = [1.0, 0.0]
@time psi = MPS(spin_up, L)
println("[+] Constructed wavefunction.")

println("[ ] Sweeping.")
outputs = SweepOutput[]
if !isnothing(output_file)
    push!(outputs, SweepOutputFile(output_file))
end
push!(outputs, SweepOutputDynamic())
@time hist = dmrg!(psi, H, SweepSchedule(max_sweeps; tolerance=5e-12), outputs=outputs)
println("[+] Swept.")

hist.converged || @warn("Ground state not converged.")
if hist[end].mps_max_bond_dim >= hist[end].set.bond_max
    @warn("Bond dimension saturated at final sweep.")
end

# Ground state energy.
println_result("E0 = $(hist[end].energy)")
# Von Neumann entanglement entropy.
println_result("SvN = $(S_vn(hist[end].middle_eigvals))")
