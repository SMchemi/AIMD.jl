# analysis.jl
module Analysis

using Plots
using Dates
using Printf
export plot_energies, write_trajectory_xyz, write_trajectory_vel, 
initialize_trajectory_file, append_trajectory_xyz, append_trajectory_vel,
initialize_md_log, append_md_log, finalize_md_log

"""
    plot_energies(times, Es_pot, Es_kin, Es_tot, filename="energies.png")

Plots potential, kinetic, and total energies vs. time and saves the plot.
"""
function plot_energies(times, Es_pot, Es_kin, Es_tot; filename="energies.png")
    if isempty(times) || isempty(Es_pot) || isempty(Es_kin) || isempty(Es_tot)
        @warn "No energy data to plot."
        return
    end
    plot(times, Es_pot, label="Potential Energy", xlabel="Time (a.u.)", ylabel="Energy (Hartree)")
    plot!(times, Es_kin, label="Kinetic Energy")
    plot!(times, Es_tot, label="Total Energy", linewidth=2)
    title!("Energy Conservation")
    savefig(filename)
    println("Energy plot saved to $filename")
end

"""
    write_trajectory_xyz(filepath::String, atom_symbols::Vector{String}, trajectory_positions::Vector{Matrix{Float64}}, times::Vector{Float64})

Writes the trajectory to an XYZ file. Converts Bohr back to Angstrom.
trajectory_positions is a Vector where each element is a [3 x N_atoms] matrix for one frame.
"""
function write_trajectory_xyz(filepath::String, atom_symbols::Vector{String}, trajectory_positions::Vector{Matrix{Float64}}, times::Vector{Float64})
    if isempty(trajectory_positions)
        @warn "No trajectory data to write."
        return
    end
    n_atoms = length(atom_symbols)
    n_frames = length(trajectory_positions)
    bohr_to_angstrom = 0.5291772109 # Bohr / Angstrom

    open(filepath, "w") do io
        for i in 1:n_frames
            println(io, n_atoms)
            @printf(io, "Frame: %d, Time: %.4f fs\n", i, times[i])
            positions_bohr = trajectory_positions[i] # [3 x N_atoms]
            for j in 1:n_atoms
                x_ang = positions_bohr[1, j] * bohr_to_angstrom
                y_ang = positions_bohr[2, j] * bohr_to_angstrom
                z_ang = positions_bohr[3, j] * bohr_to_angstrom
                @printf(io, "%-3s %15.8f %15.8f %15.8f\n", atom_symbols[j], x_ang, y_ang, z_ang)
            end
        end
    end
    println("Trajectory saved to $filepath")
end

"""
    write_trajectory_vel(filepath::String, atom_symbols::Vector{String}, trajectory_positions::Vector{Matrix{Float64}}, times::Vector{Float64})

Writes the trajectory to an XYZ file. Converts Bohr back to Angstrom.
trajectory_positions is a Vector where each element is a [3 x N_atoms] matrix for one frame.
"""
function write_trajectory_vel(filepath::String, atom_symbols::Vector{String}, trajectory_velocities::Vector{Matrix{Float64}}, times::Vector{Float64})
    if isempty(trajectory_velocities)
        @warn "No trajectory data to write."
        return
    end
    n_atoms = length(atom_symbols)
    n_frames = length(trajectory_velocities)
    bohr_to_angstrom = 0.5291772109 # Bohr / Angstrom
    fs_to_au = 41.341374575751 # fs / au 

    open(filepath, "w") do io
        for i in 1:n_frames
            println(io, n_atoms)
            @printf(io, "Frame: %d, Time: %.4f fs\n", i, times[i])
            velocities_bohr = trajectory_velocities[i] # [3 x N_atoms]
            for j in 1:n_atoms
                x_ang = velocities_bohr[1, j] * bohr_to_angstrom * fs_to_au
                y_ang = velocities_bohr[2, j] * bohr_to_angstrom * fs_to_au
                z_ang = velocities_bohr[3, j] * bohr_to_angstrom * fs_to_au
                @printf(io, "%-3s %15.8f %15.8f %15.8f\n", atom_symbols[j], x_ang, y_ang, z_ang)
            end
        end
    end
    println("Trajectory saved to $filepath")
end

"""
    initialize_trajectory_file(filepath::String, atom_symbols::Vector{String})
Creates a new trajectory file at the start of the MD simulation.
"""
function initialize_trajectory_file(filepath::String)
    # Create a new file (or overwrite existing)
    open(filepath, "w") do io
        # Just create the file, no content needed initially
    end
    println("Initialized trajectory file: $filepath")
end

"""
    append_trajectory_xyz(filepath::String, atom_symbols::Vector{String}, positions::Matrix{Float64}, step::Int, time::Float64)
Appends the current positions to an XYZ file. Converts Bohr to Angstrom.
positions is a [3 x N_atoms] matrix for the current frame.
"""
function append_trajectory_xyz(filepath::String, atom_symbols::Vector{String}, positions::Matrix{Float64}, step::Int, time::Float64)
    n_atoms = length(atom_symbols)
    bohr_to_angstrom = 0.5291772109 # Bohr / Angstrom
    
    open(filepath, "a") do io
        println(io, n_atoms)
        @printf(io, "Frame: %d, Time: %.4f fs\n", step, time)
        
        for j in 1:n_atoms
            x_ang = positions[1, j] * bohr_to_angstrom
            y_ang = positions[2, j] * bohr_to_angstrom
            z_ang = positions[3, j] * bohr_to_angstrom
            @printf(io, "%-3s %15.8f %15.8f %15.8f\n", atom_symbols[j], x_ang, y_ang, z_ang)
        end
    end
end

"""
    append_trajectory_vel(filepath::String, atom_symbols::Vector{String}, velocities::Matrix{Float64}, step::Int, time::Float64)
Appends the current velocities to a file. Converts units from atomic to Angstrom/fs.
velocities is a [3 x N_atoms] matrix for the current frame.
"""
function append_trajectory_vel(filepath::String, atom_symbols::Vector{String}, velocities::Matrix{Float64}, step::Int, time::Float64)
    n_atoms = length(atom_symbols)
    bohr_to_angstrom = 0.5291772109 # Bohr / Angstrom
    fs_to_au = 41.341374575751 # fs / au
    
    open(filepath, "a") do io
        println(io, n_atoms)
        @printf(io, "Frame: %d, Time: %.4f fs\n", step, time)
        
        for j in 1:n_atoms
            x_ang = velocities[1, j] * bohr_to_angstrom * fs_to_au
            y_ang = velocities[2, j] * bohr_to_angstrom * fs_to_au
            z_ang = velocities[3, j] * bohr_to_angstrom * fs_to_au
            @printf(io, "%-3s %15.8f %15.8f %15.8f\n", atom_symbols[j], x_ang, y_ang, z_ang)
        end
    end
end

"""
    initialize_md_log(filepath::String, config::Dict)
Creates a new MD log file with header information about the simulation setup.
"""
function initialize_md_log(filepath::String, config::Dict, atom_symbols::Vector{String}, masses::Vector{Float64})
    open(filepath, "w") do io
        println(io, "===============================================")
        println(io, "  Molecular Dynamics Simulation Log")
        println(io, "===============================================")
        println(io, "Date: ", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
        println(io, "")
        
        println(io, "Simulation Parameters:")
        println(io, "---------------------")
        println(io, "Time step: $(config[:dt_fs]) fs ($(config[:dt]) a.u.)")
        println(io, "Number of steps: $(config[:n_steps])")
        println(io, "Target state: $(config[:target_state_index])")
        println(io, "Save frequency: Every $(config[:save_freq]) steps")
        println(io, "")
        
        println(io, "System Information:")
        println(io, "------------------")
        println(io, "Number of atoms: $(length(atom_symbols))")
        println(io, "Atoms: $(join(atom_symbols, ", "))")
        println(io, "")
        
        println(io, "Initial Files:")
        println(io, "-------------")
        println(io, "Position file: $(config[:position_file])")
        println(io, "Velocity file: $(config[:velocity_file])")
        println(io, "")
        
        println(io, "===============================================")
        println(io, "            Simulation Progress")
        println(io, "===============================================")
        println(io, "Step      Time(fs)    E_pot(H)      E_kin(H)      E_tot(H)      ΔE_tot      Temperature(K)")
        println(io, "-------------------------------------------------------------------------------------")
    end
    println("MD log file initialized: $filepath")
end

"""
    append_md_log(filepath::String, step::Int, time::Float64, e_pot::Float64, e_kin::Float64, e_tot::Float64, 
                  initial_etot::Float64, masses::Vector{Float64}, velocities::Matrix{Float64}, notes::String="")
Appends a line to the MD log file with current simulation data.
Optional 'notes' parameter allows adding specific events or warnings.
"""
function append_md_log(filepath::String, step::Int, time::Float64, e_pot::Float64, e_kin::Float64, e_tot::Float64, 
                       initial_etot::Float64, masses::Vector{Float64}, velocities::Matrix{Float64}, notes::String="")
    
    # Calculate temperature from kinetic energy (if needed)
    # T = 2*E_kin / (3*N*k_B), where k_B is Boltzmann constant
    # In atomic units: k_B = 3.1668114e-6 Ha/K
    k_B_au = 3.1668114e-6
    n_atoms = length(masses)
    temperature = 2.0 * e_kin / (3.0 * n_atoms * k_B_au)
    
    # Calculate energy drift as percentage of initial total energy
    energy_drift = e_tot - initial_etot
    energy_drift_pct = 100.0 * energy_drift / abs(initial_etot)
    
    open(filepath, "a") do io
        @printf(io, "%-8d  %10.4f  %12.8f  %12.8f  %12.8f  %+8.2e  %10.2f", 
                step, time, e_pot, e_kin, e_tot, energy_drift, temperature)
        
        # Add any notes
        if !isempty(notes)
            print(io, "  # $notes")
        end
        println(io)
        
        # Add extra details for milestone steps
        if step % 1000 == 0
            println(io, "  -- 1000-step milestone reached --")
        end
    end
end

"""
    finalize_md_log(filepath::String, total_steps::Int, elapsed_time::Float64, final_etot::Float64, initial_etot::Float64)
Writes final summary information to the MD log file.
"""
function finalize_md_log(filepath::String, total_steps::Int, elapsed_time::Float64, final_etot::Float64, initial_etot::Float64)
    energy_drift = final_etot - initial_etot
    energy_drift_pct = 100.0 * energy_drift / abs(initial_etot)
    
    open(filepath, "a") do io
        println(io, "")
        println(io, "===============================================")
        println(io, "            Simulation Summary")
        println(io, "===============================================")
        println(io, "Total steps completed: $total_steps")
        println(io, "Wall clock time: $elapsed_time seconds")
        @printf(io, "Initial total energy: %.8f Hartree\n", initial_etot)
        @printf(io, "Final total energy: %.8f Hartree\n", final_etot)
        @printf(io, "Energy drift: %.8f Hartree (%.6f%%)\n", energy_drift, energy_drift_pct)
        
        # Energy conservation assessment
        if abs(energy_drift_pct) < 0.01
            println(io, "Energy conservation: EXCELLENT")
        elseif abs(energy_drift_pct) < 0.1
            println(io, "Energy conservation: GOOD")
        elseif abs(energy_drift_pct) < 1.0
            println(io, "Energy conservation: ACCEPTABLE")
        else
            println(io, "Energy conservation: POOR - Check time step or integration method")
        end
        
        println(io, "")
        println(io, "Simulation completed at: ", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
        println(io, "===============================================")
    end
    println("MD log finalized: $filepath")
end

end # module Analysis
