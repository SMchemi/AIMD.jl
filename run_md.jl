# run_md.jl
using LinearAlgebra # For matrix operations
using Printf      # For formatted output
using Plots       # For plotting energies
# Choose a backend for Plots, e.g., GR. You might need to add Plots and GR packages.
# using Pkg; Pkg.add("Plots"); Pkg.add("GR")
gr()

# Include custom modules (assuming they are in the same directory or accessible via LOAD_PATH)
include("config.jl")
include("file_parsers.jl")
include("orca_interface.jl")
include("md_integrators.jl")
include("analysis.jl")
include("wigner_sampling.jl")

# Use functions from modules
using .SimulationConfig
using .FileParsers
using .WignerSampling
using .OrcaInterface
using .MDIntegrators
using .Analysis

original_dir = "/home/smchemi/julia/MD_Interface"

function run_simulation()
    # --- 1. Initialization & Setup ---
    println("--- Initializing MD Simulation ---")
    cfg = SimulationConfig.settings # Get config dictionary

    # --- Create Base Output Directory and Change Into It ---
    base_output_dir = abspath(original_dir, cfg[:output_dir])
    if !isdir(base_output_dir)
        println("Creating base output directory: $base_output_dir")
        mkpath(base_output_dir)
    end
    # Change the current working directory to the base output directory
    # This simplifies relative path handling for step directories and final output
    try
        cd(base_output_dir)
        println("Changed working directory to: $(pwd())")
    catch e
        error("Could not change directory to output directory '$base_output_dir': $e")
    end
    # ------------------------------------------------------

    # --- Load Initial Conditions ---
    # Paths in config are now interpreted relative to the original launch directory
    # OR should be absolute paths. We use abspath to handle both cases if config paths are relative.
    initial_pos_file_path = abspath(original_dir, cfg[:position_file]) # Use original_dir if needed, or just cfg[:position_file] if absolute
    initial_vel_file_path = abspath(original_dir, cfg[:velocity_file]) # Assuming original_dir holds the path where Julia was launched

    println("Reading initial positions from: $initial_pos_file_path")
    atom_symbols, initial_positions = FileParsers.read_xyz(initial_pos_file_path)
    n_atoms = length(atom_symbols)
    println("Number of atoms: $n_atoms")

    if length(cfg[:masses]) != n_atoms
        error("Number of masses ($(length(cfg[:masses]))) does not match number of atoms ($n_atoms)")
    end
    # Masses should be in atomic units (electron mass = 1) as per config.jl
    masses_au = cfg[:masses]

    println("Reading initial velocities from: $initial_vel_file_path")
    initial_velocities = FileParsers.read_velocities(initial_vel_file_path, n_atoms)
    println("Initial velocities loaded successfully.")
    # Positions/Velocities are typically read in Angstrom or Angstrom/fs.
    # MD uses Bohr and atomic units. Conversion needed? Assuming units match MD integrator for now.
    # Let's assume MDIntegrators expect Bohr positions and Bohr/a.u. time velocities.
    # TODO: Add unit conversions here if initial files are not in atomic units.
    initial_positions_bohr = initial_positions ./ 0.529177210903
    initial_velocities_au = initial_velocities ./ 0.529177210903 ./ 41.341374575751 # (angstrom_to_bohr/ fs_to_au_time)

    # --- Simulation Parameters ---
    dt = cfg[:dt]            # Time step in atomic units of time
    n_steps = cfg[:n_steps]
    target_state = cfg[:target_state_index]

    # --- Initialize State Variables ---
    # Use atomic units internally for MD
    # Assuming initial files were converted if necessary
    current_positions = copy(initial_positions_bohr) # Assume Bohr
    current_velocities = copy(initial_velocities_au) # Assume Bohr / a.u. time
    current_forces = zeros(Float64, 3, n_atoms) # Hartree / Bohr (from ORCA)
    current_potential_energy = 0.0 # Hartree
    current_kinetic_energy = 0.0   # Hartree
    current_total_energy = 0.0     # Hartree
    time = 0.0                     # Atomic units of time

    # --- Setup Storage Arrays ---
    save_freq = cfg[:save_freq]
    num_saves = floor(Int, n_steps / save_freq) + 1 # Number of frames including initial
    save_indices = [0; [i * save_freq for i in 1:floor(Int, n_steps / save_freq)]] # Steps 0, 10, 20...
    if n_steps % save_freq != 0 # Ensure last step is saved if not on interval
         push!(save_indices, n_steps)
         save_indices = unique(sort(save_indices))
         num_saves = length(save_indices)
    else
        num_saves = length(save_indices) # Already includes 0
    end

    times_stored = Vector{Float64}(undef, num_saves)
    # Store positions in Bohr for trajectory file consistency (can convert later if needed)
    # trajectory_positions_bohr = Vector{Matrix{Float64}}(undef, num_saves)
    # trajectory_velocities_au = Vector{Matrix{Float64}}(undef, num_saves)
    potential_energies = Vector{Float64}(undef, num_saves)
    kinetic_energies = Vector{Float64}(undef, num_saves)
    total_energies = Vector{Float64}(undef, num_saves)
    save_counter = 0 # Index for storage arrays

    # --- 2. Initial Force Calculation ---
    println("--- Calculating Initial Forces (Step 0) using Direct Execution / .engrad ---")
    # Pass positions in Bohr, get forces in Hartree/Bohr
    current_potential_energy, current_forces = OrcaInterface.calculate_orca_energy_forces_engrad(
        current_positions, atom_symbols, cfg, 0 # Step 0
    )
    current_kinetic_energy = MDIntegrators.calculate_kinetic_energy(masses_au, current_velocities)
    current_total_energy = current_potential_energy + current_kinetic_energy

    # Store initial state (Step 0) - check if 0 is in save_indices (it should be)
    if 0 in save_indices
        save_counter += 1
        times_stored[save_counter] = time
        # trajectory_positions_bohr[save_counter] = copy(current_positions)
	    # trajectory_velocities_au[save_counter] = copy(current_velocities)
        potential_energies[save_counter] = current_potential_energy
        kinetic_energies[save_counter] = current_kinetic_energy
        total_energies[save_counter] = current_total_energy
    else
         @warn "Step 0 was not added to save_indices. Check save logic."
    end

    println("Initial State (Time: $time a.u.):")
    @printf("  E_pot: %.8f Hartree\n", current_potential_energy)
    @printf("  E_kin: %.8f Hartree\n", current_kinetic_energy)
    @printf("  E_tot: %.8f Hartree\n", current_total_energy)
    println("--- Starting MD Loop ---")

    # Define output file paths (relative to the current directory, which is base_output_dir)
    energy_plot_file = "$(cfg[:output_prefix])_energies.png"
    trajectory_file = "$(cfg[:output_prefix])_trajectory.xyz"
    velocity_file = "$(cfg[:output_prefix])_velocity.dat"
    log_file = "$(cfg[:output_prefix])_md.log"
    start_time = time_ns() # For tracking wall clock time
    Analysis.initialize_md_log(log_file, cfg, atom_symbols, masses_au)
    
    # Log initial state
    initial_total_energy = current_total_energy
    Analysis.append_md_log(log_file, 0, 0.0, current_potential_energy, current_kinetic_energy, 
                          current_total_energy, initial_total_energy, masses_au, current_velocities, 
                          "Initial configuration")
    Analysis.initialize_trajectory_file(trajectory_file)
    Analysis.initialize_trajectory_file(velocity_file)

    # --- 3. Main MD Loop ---
    for step = 1:n_steps
        # --- Velocity Verlet Step ---
        # 1. Velocity Update (Half Step) & Position Update using previous forces F(t)
        #    Updates positions to r(t+dt) and velocities to v(t+dt/2)
        #    Requires inputs in atomic units: Bohr, Bohr/(a.u. time), Hartree/Bohr, a.u. mass, a.u. time
        MDIntegrators.velocity_verlet_step!(current_positions, current_velocities, current_forces, masses_au, dt)

        # 3. Force/Energy Update (Call ORCA for the new positions r(t+dt))
        #    Calculates F(t+dt) and E_pot(t+dt)
        new_potential_energy, new_forces = OrcaInterface.calculate_orca_energy_forces_engrad(
            current_positions, atom_symbols, cfg, step
        )

        # 4. Velocity Update (Second Half Step) using *new* forces F(t+dt)
        #    Calculates v(t+dt) = v(t+dt/2) + 0.5 * a(t+dt) * dt
        #    Requires forces in Hartree/Bohr, mass in a.u. mass, dt in a.u. time
        mass_reshaped = reshape(masses_au, 1, n_atoms)
        new_acceleration = new_forces ./ mass_reshaped # Bohr / (a.u. time)^2
        current_velocities .+= 0.5 .* new_acceleration .* dt # Update v to v(t+dt)
        # --- End Velocity Verlet Step ---

        # 5. Update State for next iteration
        current_forces = new_forces
        current_potential_energy = new_potential_energy

        # 6. Calculate Kinetic Energy at t+dt
        current_kinetic_energy = MDIntegrators.calculate_kinetic_energy(masses_au, current_velocities)

        # 7. Calculate Total Energy at t+dt
        current_total_energy = current_potential_energy + current_kinetic_energy

        # 8. Update Time
        time += dt

        # 9. Store Data (if current step is a save step)
        if step in save_indices
            # Write to files at the specified save frequency
            Analysis.append_trajectory_xyz(trajectory_file, atom_symbols, current_positions, step, time*cfg[:dt_fs])
            Analysis.append_trajectory_vel(velocity_file, atom_symbols, current_velocities, step, time*cfg[:dt_fs])
    
            save_counter += 1
            if save_counter <= num_saves
                times_stored[save_counter] = time
                # trajectory_positions_bohr[save_counter] = copy(current_positions)
		        # trajectory_velocities_au[save_counter] = copy(current_velocities)
                potential_energies[save_counter] = current_potential_energy
                kinetic_energies[save_counter] = current_kinetic_energy
                total_energies[save_counter] = current_total_energy
            else
                @warn "Save counter ($save_counter) exceeds allocated size ($num_saves) at step $step. Check save logic."
            end
        end

        # 10. Monitor Progress
        if step % cfg[:print_freq] == 0 || step == n_steps
            # Print to console
            @printf("Step: %d / %d | Time: %.2f fs | E_pot: %.8f | E_kin: %.8f | E_tot: %.8f\n",
                    step, n_steps, step*cfg[:dt_fs], current_potential_energy, current_kinetic_energy, current_total_energy)
            
            # Log to file
            Analysis.append_md_log(log_file, step, step*cfg[:dt_fs], current_potential_energy, current_kinetic_energy, 
                                current_total_energy, initial_total_energy, masses_au, current_velocities)
        end

    end # End MD Loop

    println("--- MD Loop Finished ---")

    # Trim storage arrays if the last step wasn't saved (shouldn't happen with current logic, but safe check)
    if save_counter < num_saves
        @warn "Trimming storage arrays, expected $num_saves saves, got $save_counter."
        times_stored = times_stored[1:save_counter]
        # trajectory_positions_bohr = trajectory_positions_bohr[1:save_counter]
	    # trajectory_velocities_au = trajectory_velocities_bohr[1:save_counter]
        potential_energies = potential_energies[1:save_counter]
        kinetic_energies = kinetic_energies[1:save_counter]
        total_energies = total_energies[1:save_counter]
    end

    # --- 4. Analysis & Finalization ---
    println("--- Post-processing and Analysis ---")

    elapsed_time = (time_ns() - start_time) / 1.0e9  # Convert ns to seconds
    Analysis.finalize_md_log(log_file, n_steps, elapsed_time, current_total_energy, initial_total_energy)

    # Plot Energies (times are in a.u., energies in Hartree)
    Analysis.plot_energies(collect(0:n_steps)*cfg[:dt_fs], potential_energies, kinetic_energies, total_energies, filename=energy_plot_file)

    # Write Trajectory (positions are stored in Bohr)
    # write_trajectory_xyz likely expects Angstrom. Convert Bohr positions to Angstrom before writing.
    # Analysis.write_trajectory_xyz(trajectory_file, atom_symbols, trajectory_positions_bohr, collect(0:n_steps)*cfg[:dt_fs])
    # Analysis.write_trajectory_vel(velocity_file, atom_symbols, trajectory_velocities_au, collect(0:n_steps)*cfg[:dt_fs])
    # ^^^ NOTE: Added convert_to_angstrom flag assumption to analysis function

    println("--- Simulation Complete ---")
    println("Output files saved in: $(pwd())") # pwd() is base_output_dir

end # function run_simulation

# --- Execute the Simulation ---
run_simulation()
