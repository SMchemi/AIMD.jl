# wigner_sampling.jl
# Implementation of Wigner sampling for molecular dynamics initial conditions

module WignerSampling

export generate_wigner_samples, sample_wigner_from_hess, read_vibrational_data, read_xyz_structure, write_md_initial_conditions

using LinearAlgebra
using Printf

"""
    read_vibrational_data(hess_filename::String)

Read vibrational frequencies and normal modes from an ORCA .hess file.
Returns frequencies, normal modes, masses, and number of atoms.
"""
function read_vibrational_data(hess_filename::String)
    open(hess_filename, "r") do f
        lines = readlines(f)
        
        # Find vibrational frequencies
        freq_start = findfirst(l -> occursin("\$vibrational_frequencies", l), lines) + 1
        n_freq = parse(Int, lines[freq_start])
        freqs = Float64[]
        for i in 1:n_freq
            freq = parse(Float64, split(lines[freq_start+i])[1])
            push!(freqs, freq)
        end
        
        # Find normal modes
        modes_start = findfirst(l -> occursin("\$normal_modes", l), lines) + 1
        dims = split(lines[modes_start])
        n_modes = parse(Int, dims[1])
        
        # Read normal mode matrix
        normal_modes = zeros(Float64, n_modes, n_modes)
        mode_idx = 0
        line_idx = modes_start + 1
        
        while mode_idx < n_modes
            # Skip header
            line_idx += 1
            
            for i in 1:n_modes
                values = split(lines[line_idx+i])
                for j in 1:min(5, n_modes-mode_idx)
                    normal_modes[i, mode_idx+j] = parse(Float64, values[j])
                end
            end
            
            mode_idx += 5
            line_idx += n_modes + 1
        end
        
        # Read atomic masses
        atoms_start = findfirst(l -> occursin("\$atoms", l), lines) + 1
        n_atoms = parse(Int, lines[atoms_start])
        masses = Float64[]
        atom_symbols = String[]
        
        for i in 1:n_atoms
            parts = split(lines[atoms_start+i])
            push!(atom_symbols, parts[1])
            mass = parse(Float64, parts[2])
            push!(masses, mass)
        end
        
        return freqs, normal_modes, masses, atom_symbols, n_atoms
    end
end

"""
    read_xyz_structure(xyz_filename::String)

Read molecular structure from an XYZ file.
Returns atomic symbols and positions (in Angstroms).
"""
function read_xyz_structure(xyz_filename::String)
    open(xyz_filename, "r") do f
        lines = readlines(f)
        
        # First line contains the number of atoms
        n_atoms = parse(Int, lines[1])
        
        # Second line is usually a comment
        positions = zeros(Float64, 3, n_atoms)
        atom_symbols = String[]
        
        # From the third line, atom coordinates
        for i in 1:n_atoms
            parts = split(lines[i+2])
            push!(atom_symbols, parts[1])
            positions[1, i] = parse(Float64, parts[2])
            positions[2, i] = parse(Float64, parts[3])
            positions[3, i] = parse(Float64, parts[4])
        end
        
        return atom_symbols, positions
    end
end

"""
    coth(x::Float64)

Calculate the hyperbolic cotangent function.
"""
function coth(x::Float64)
    if abs(x) < 1e-10
        return 1.0/x
    else
        return (exp(2x) + 1) / (exp(2x) - 1)
    end
end

"""
    generate_wigner_samples(freqs::Vector{Float64}, normal_modes::Matrix{Float64}, 
                           masses::Vector{Float64}, eq_positions::Matrix{Float64}, 
                           T_K::Float64, n_samples::Int)

Generate Wigner distribution samples for molecular dynamics initial conditions.
Returns positions (in Bohr) and velocities (in atomic units) for each sample.
"""
function generate_wigner_samples(freqs::Vector{Float64}, normal_modes::Matrix{Float64}, 
                                masses::Vector{Float64}, eq_positions::Matrix{Float64}, 
                                T_K::Float64, n_samples::Int)
    # Constants
    kb = 3.166811563e-6  # Hartree/K (Boltzmann constant)
    hbar = 1.0  # ħ = 1 in atomic units
    
    # 1. Exclude the 6 modes with smallest absolute values
    abs_freqs = abs.(freqs)
    sorted_idx = sortperm(abs_freqs)
    excluded_idx = sorted_idx[1:6]  # Smallest 6 modes
    
    # 2. From remaining modes, exclude imaginary frequencies (negative values)
    remaining_idx = setdiff(1:length(freqs), excluded_idx)
    imag_idx = filter(i -> freqs[i] < 0, remaining_idx)
    
    # 3. Final set of modes to use (remaining positive frequencies)
    vib_idx = setdiff(remaining_idx, imag_idx)
    
    # Storage arrays
    n_atoms = length(masses)
    # Note: Positions are stored as [xyz, atom, sample]
    positions = zeros(Float64, 3, n_atoms, n_samples)
    momenta = zeros(Float64, 3, n_atoms, n_samples)
    
    # Set equilibrium positions for all samples
    for s in 1:n_samples
        positions[:, :, s] = eq_positions
    end
    
    # Sample from Wigner distribution for each mode
    for s in 1:n_samples
        for idx in vib_idx
            freq = freqs[idx]
            
            # Convert frequency to angular frequency (cm^-1 to atomic units)
            ω = freq * 4.556335e-6  # 2π·c·ν in atomic units
            
            # Calculate thermal energy factor
            beta = 1.0 / (kb * T_K)
            coth_factor = coth(beta * hbar * ω / 2)
            
            # Standard deviations for Gaussian sampling 
            σ_q = sqrt(1 / (2 * ω))
            σ_p = sqrt(ω / 2)
            
            # Adjust for temperature effects
            σ_q *= sqrt(coth_factor)
            σ_p *= sqrt(coth_factor)
            
            # Sample from normal distribution
            q_sample = randn() * σ_q
            p_sample = randn() * σ_p
            
            # Transform normal mode coordinates to Cartesian coordinates
            for a in 1:n_atoms
                mass_sqrt = sqrt(masses[a])
                for xyz in 1:3
                    dof = (a-1)*3 + xyz
                    mode_coeff = normal_modes[dof, idx]
                    
                    # Apply mass weighting
                    positions[xyz, a, s] += q_sample * mode_coeff / mass_sqrt
                    momenta[xyz, a, s] += p_sample * mode_coeff * mass_sqrt
                end
            end
        end
    end
    
    # Convert momenta to velocities
    velocities = similar(momenta)
    for s in 1:n_samples
        for a in 1:n_atoms
            velocities[:, a, s] = momenta[:, a, s] / masses[a]
        end
    end
    
    return positions, velocities
end

"""
    sample_wigner_from_hess(hess_file::String, xyz_file::String, T_K::Float64, n_samples::Int)

Generate Wigner samples using frequency data from a .hess file and structure from an .xyz file.
Returns atomic symbols, positions (in Bohr), velocities (in atomic units), and masses.
"""
function sample_wigner_from_hess(hess_file::String, xyz_file::String, T_K::Float64, n_samples::Int)
    # Read vibrational data
    freqs, normal_modes, masses, atom_symbols_hess, n_atoms = read_vibrational_data(hess_file)
    
    # Read equilibrium geometry
    atom_symbols_xyz, eq_positions = read_xyz_structure(xyz_file)
    
    # Ensure atom counts match
    if length(atom_symbols_xyz) != n_atoms
        error("Number of atoms in XYZ file ($(length(atom_symbols_xyz))) doesn't match HESS file ($n_atoms)")
    end
    
    # Convert positions from Angstrom to Bohr
    eq_positions_bohr = eq_positions ./ 0.529177210903
    
    # Generate Wigner samples
    positions_bohr, velocities_au = generate_wigner_samples(
        freqs, normal_modes, masses, eq_positions_bohr, T_K, n_samples
    )
    
    return atom_symbols_hess, positions_bohr, velocities_au, masses
end

"""
    write_md_initial_conditions(atom_symbols::Vector{String}, positions::Array{Float64,3}, 
                               velocities::Array{Float64,3}, sample_idx::Int, 
                               positions_file::String, velocities_file::String)

Write the initial positions and velocities for MD simulation to files.
Positions are written in Angstrom, velocities in Angstrom/fs.
"""
function write_md_initial_conditions(atom_symbols::Vector{String}, positions::Array{Float64,3}, 
                                    velocities::Array{Float64,3}, sample_idx::Int, 
                                    positions_file::String, velocities_file::String)
    n_atoms = length(atom_symbols)
    
    # Extract the selected sample
    selected_positions = positions[:, :, sample_idx]  # [3, n_atoms]
    selected_velocities = velocities[:, :, sample_idx]  # [3, n_atoms]
    
    # Convert from atomic units to Angstrom and Angstrom/fs
    bohr_to_angstrom = 0.529177210903
    au_time_to_fs = 0.02418884254  # 1 a.u. of time = 0.02419 fs
    
    positions_angstrom = selected_positions .* bohr_to_angstrom
    velocities_angstrom_fs = selected_velocities .* bohr_to_angstrom ./ au_time_to_fs
    
    # Write positions to XYZ file
    open(positions_file, "w") do f
        println(f, n_atoms)
        println(f, "Initial positions for MD - Wigner sample $sample_idx")
        for i in 1:n_atoms
            @printf(f, "%-2s %14.8f %14.8f %14.8f\n", 
                   atom_symbols[i], positions_angstrom[1, i], positions_angstrom[2, i], positions_angstrom[3, i])
        end
    end
    
    # Write velocities to velocity file
    open(velocities_file, "w") do f
        println(f, n_atoms)
        println(f, "Initial velocities for MD - Wigner sample $sample_idx (Angstrom/fs)")
        for i in 1:n_atoms
            @printf(f, "%-2s %14.8f %14.8f %14.8f\n", 
                   atom_symbols[i], velocities_angstrom_fs[1, i], velocities_angstrom_fs[2, i], velocities_angstrom_fs[3, i])
        end
    end
    
    println("Initial conditions written to:")
    println("  Positions: $positions_file")
    println("  Velocities: $velocities_file")
end

end  # module WignerSampling