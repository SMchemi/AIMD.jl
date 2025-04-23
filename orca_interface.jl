# orca_interface.jl
module OrcaInterface

using ..SimulationConfig
using Printf
using LinearAlgebra

export calculate_orca_energy_forces_engrad # Export the new function name

# --- Keep coordinate formatter ---
"""
Formats coordinates (Bohr) into ORCA's XYZ format string (Angstrom).
"""
function format_coords_for_orca(atom_symbols::Vector{String}, positions_bohr::Matrix{Float64})
    n_atoms = length(atom_symbols)
    bohr_to_angstrom = 0.529177210903
    lines = []
    for i in 1:n_atoms
        pos_ang = positions_bohr[:, i] .* bohr_to_angstrom
        push!(lines, @sprintf("%-3s %16.10f %16.10f %16.10f",
                             atom_symbols[i], pos_ang[1], pos_ang[2], pos_ang[3]))
    end
    return join(lines, "\n")
end


"""
Parses ORCA .engrad file for energy and gradient.
Returns: potential_energy (Float64, Hartree), gradient (Matrix{Float64}, 3xN, Hartree/Bohr)
"""
function parse_orca_engrad(engrad_file::String)
    if !isfile(engrad_file)
        @warn "Cannot parse .engrad file: File not found at $engrad_file"
        return NaN, nothing
    end
    
    lines = readlines(engrad_file)
    n_atoms = -1
    energy = NaN
    gradient = nothing
    
    try
        # Find number of atoms
        idx_natoms_header = findfirst(l -> occursin("# Number of atoms", l), lines)
        if idx_natoms_header === nothing || idx_natoms_header + 1 > length(lines)
            error("Could not find '# Number of atoms' header or value.")
        end
        n_atoms = parse(Int, strip(lines[idx_natoms_header + 2]))
        if n_atoms <= 0 
            error("Parsed invalid number of atoms: $n_atoms") 
        end
        
        # Find energy
        idx_energy_header = findfirst(l -> occursin("# The current total energy in Eh", l), lines)
        if idx_energy_header === nothing || idx_energy_header + 1 > length(lines)
            error("Could not find '# The current total energy in Eh' header or value.")
        end
        energy = parse(Float64, strip(lines[idx_energy_header + 2]))
        
        # Find gradient block
        idx_grad_header = findfirst(l -> occursin("# The current gradient in Eh/bohr", l), lines)
        if idx_grad_header === nothing
            error("Could not find '# The current gradient in Eh/bohr' header.")
        end
        
        grad_start_line = idx_grad_header + 2
        num_grad_components = 3 * n_atoms
        
        # Filter out any lines that don't contain numerical data
        potential_grad_lines = lines[grad_start_line:end]
        gradient_lines = String[]
        
        for line in potential_grad_lines
            stripped = strip(line)
            # Skip lines that are empty or contain a comment marker
            if isempty(stripped) || startswith(stripped, "#")
                continue
            end
            
            # Try to parse as float to verify it's a gradient line
            try
                parse(Float64, stripped)
                push!(gradient_lines, stripped)
                # Break if we've collected enough gradient components
                if length(gradient_lines) >= num_grad_components
                    break
                end
            catch
                # Not a valid floating-point number, might be the start of a new section
                break
            end
        end
        
        # Check if we found the right number of gradient components
        if length(gradient_lines) != num_grad_components
            error("Expected $num_grad_components gradient components, found $(length(gradient_lines)).")
        end
        
        # Parse gradient lines and reshape
        gradient_vector = parse.(Float64, gradient_lines)
        gradient = reshape(gradient_vector, 3, n_atoms) # Reshape into 3xN matrix
    catch e
        @warn "Error parsing .engrad file '$engrad_file': $e"
        return NaN, nothing # Return failure indication
    end
    
    # Final checks
    if isnan(energy) || gradient === nothing || size(gradient) != (3, n_atoms)
        @warn "Parsing .engrad file '$engrad_file' resulted in incomplete data."
        return NaN, nothing
    end
    
    return energy, gradient # Gradient G (Hartree/Bohr)
end

"""
Runs ORCA directly within a step-specific subdirectory of the main output folder.
Parses the .engrad file for results. Manages .gbw continuity.
"""
function calculate_orca_energy_forces_engrad(positions_bohr::Matrix{Float64}, atom_symbols::Vector{String}, cfg::Dict, step::Int)
    original_dir = "/home/smchemi/julia/MD_Interface"

    try
        # --- Get Configuration ---
        template_file = abspath(original_dir, cfg[:orca_template_file])
        base_out_dir = abspath(original_dir, cfg[:output_dir]) # Base directory where run_md.jl is, contains step dirs
        target_state = cfg[:target_state_index] # Needed for warnings, not used in parsing .engrad directly
        orca_executable = cfg[:orca_executable]
        step_dir_prefix = cfg[:step_dir_prefix]
        cleanup = cfg[:cleanup_step_files]

        # --- Prepare Paths for this Step ---
	step_dir_path = "$(step_dir_prefix)_$(string(step,pad=3))"
        mkpath(step_dir_path) # Ensure step directory exists
	println(step_dir_path)

	base_filename = "orca_calc_step_$(string(step,pad=3))" # Consistent base name for files within step dir
        inp_file = joinpath(base_out_dir, step_dir_path, base_filename * ".inp")
        engrad_file = joinpath(base_out_dir, step_dir_path, base_filename * ".engrad")
        gbw_file = joinpath(base_out_dir, step_dir_path, base_filename * ".gbw") # Current step's GBW file
        out_file = joinpath(base_out_dir, step_dir_path, base_filename * ".out") # ORCA stdout log
        err_file = joinpath(base_out_dir, step_dir_path, base_filename * ".err") # ORCA stderr log

        # --- Handle Wavefunction Guess (.gbw) ---
        if step > 0
            prev_step_base = "$(step_dir_prefix)_$(string(step-1,pad=3))"
            prev_gbw_file = joinpath(base_out_dir, prev_step_base, "$(prev_step_base).gbw")
            if isfile(prev_gbw_file)
                try
                    cp(prev_gbw_file, gbw_file, force=true) # Copy previous GBW to current step dir with current name
                    println("Copied previous GBW file: $prev_gbw_file -> $gbw_file")
                catch e_cp
                    @warn "Could not copy previous GBW file from $prev_gbw_file: $e_cp"
                    # Continue without guess file, ORCA will generate anew
                end
            else
                println("No previous GBW file found at $prev_gbw_file. ORCA will generate guess.")
            end
        end

        # --- Prepare ORCA Input File Content ---
        template_content = ""
        try
            template_content = read(template_file, String)
        catch e
            error("Cannot read ORCA template file '$template_file': $e")
        end

        # Keyword checks (Engrad, iroot) - keep warnings
        if !occursin(r"!(?:.*\s)?Engrad(?:\s.*)?", template_content) && !occursin(r"!(?:.*\s)?Opt(?:\s.*)?", template_content)
            @warn "ORCA template $template_file missing '! Engrad' keyword."
        end
        if target_state > 0 && !occursin(Regex("\\biroot\\s+$target_state\\b", "i"), template_content)
            @warn "Target state is $target_state > 0, but 'iroot $target_state' not found in template $template_file."
        end

        # Format coordinates and create final input
        coords_str = format_coords_for_orca(atom_symbols, positions_bohr)
        final_inp_content = replace(template_content, "%coords%" => coords_str) # Assumes placeholder exists
        write(inp_file, final_inp_content)

        # --- Execute ORCA ---
        println("Running ORCA for step $step in $step_dir_path...")
        orca_cmd = Cmd([orca_executable, basename(inp_file)]) # Command using relative path inside step dir

        process_success = false
        run_error = nothing
        try
            cd(step_dir_path) # <<< Change directory INTO the step directory

            # Run ORCA, redirecting stdout/stderr to files *within this step directory*
            run(pipeline(orca_cmd, stdout=basename(out_file), stderr=basename(err_file)))
            process_success = true # If run completes without error

        catch e
            run_error = e
            @warn "ORCA execution failed for step $step. Error: $e"
            # process_success remains false
        finally
	    cd(base_out_dir) # <<< Change back to original directory
        end

        # --- Check Success ---
        if !process_success
            println("ORCA command failed to launch or threw error. Check files in $step_dir_path, especially $(basename(err_file))")
            throw(run_error) # Rethrow the execution error
        end

        # Primary check: Did the .engrad file get created?
        if !isfile(engrad_file)
            println("ORCA ran but the gradient file '$engrad_file' was not created.")
            println("Check ORCA output log '$out_file' and error log '$err_file' in $step_dir_path for ORCA errors.")
            # Check for common ORCA error message in .out file as a hint
            try
                if isfile(out_file)
                    out_content = read(out_file, String)
                    if occursin("INPUT ERROR", out_content) || occursin("Error in input", out_content) || occursin("TERMINATED WITHOUT FINISHING", out_content)
                        @warn "Found potential ORCA error messages in $out_file."
                    end
                end
            catch end
            error("ORCA failed to produce gradient file for step $step.")
        end

        # Optional secondary check: ORCA normal termination in .out file
        orca_terminated_normally = false
        try
            if isfile(out_file)
                open(out_file) do f
                    seekend(f); filesize = position(f); seek(f, max(0, filesize - 1024))
                    if occursin("ORCA TERMINATED NORMALLY", read(f, String))
                        orca_terminated_normally = true
                    end
                end
            end
        catch e_read
            @warn "Could not read/check ORCA output file $out_file: $e_read"
        end
        if !orca_terminated_normally
            @warn "ORCA .engrad file exists, but 'ORCA TERMINATED NORMALLY' not found in $out_file. Proceeding with parsing gradient, but check logs."
        end


        # --- Parse Results ---
        println("Parsing ORCA gradient file: $engrad_file")
        potential_energy, gradient = parse_orca_engrad(engrad_file)

        if isnan(potential_energy) || gradient === nothing
            error("Successfully ran ORCA, but failed to parse results from $engrad_file.")
        end

        forces = -gradient # F = -G (Units: Hartree/Bohr)

        # --- Cleanup Extra Files (Optional) ---
        if cleanup
            println("Cleaning up extra files in $step_dir_path...")
            files_to_keep = [basename(inp_file), basename(engrad_file), basename(gbw_file)] # Keep essential I/O
            files_to_keep_patterns = [r"\.out$", r"\.err$"] # Also keep logs by default maybe?
            try
                cd(step_dir_path)
                all_files = readdir()
                for f in all_files
                    keep = false
                    if f in files_to_keep
                        keep = true
                    else
                        for pattern in files_to_keep_patterns
                            if occursin(pattern, f)
                                keep = true
                                break
                            end
                        end
                    end

                    if !keep && isfile(f) # Delete if not marked to keep and is a file
                        try; rm(f); catch rm_err; @warn "Cleanup failed for $f: $rm_err"; end
                    end
                end
            catch e_list
                @warn "Could not list/cleanup files in $step_dir_path: $e_list"
            finally
                cd(base_out_dir) # Change back after attempting cleanup
            end
        end

        @printf("Step %d: E_pot = %.8f Hartree (from .engrad)\n", step, potential_energy)
        return potential_energy, forces

    catch e # Catch errors from the overall process
        println("\n--- ERROR during ORCA interface for step $step ---")
        println(e)
        println("-------------------------------------------------")
        # Ensure we are back in the original directory
        if pwd() != original_dir
            @warn "Attempting to return to original directory $original_dir after error."
            try; cd(original_dir); catch; end
        end
         # No complex scratch dir to clean, step dir remains for debugging unless cleanup=true ran
        rethrow(e)
    end
end


end # module OrcaInterface
