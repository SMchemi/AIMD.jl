# config.jl
module SimulationConfig

export settings

# --- User Defined Settings ---

settings = Dict(
    # File Paths
    :position_file => "initial_conditions/h2o12-0K.xyz",
    :velocity_file => "initial_conditions/h2o12-0K.dat",
    :orca_template_file => "templates/core_ext_template-14.inp",
    :output_dir => "core_14_long_output", # Base directory for MD run & step subdirs

    # MD Parameters
    :dt_fs => 0.2,
    :dt => 41.341374575751,        # Time step in atomic units of time (~0.5 fs)
    :n_steps => 100,
    :thermostat => :none,

    # System & QM Parameters
    :masses => [15.99491 * 1822.888486, 1.00784 * 1822.888486, 1.00784 * 1822.888486,
	        15.99491 * 1822.888486, 1.00784 * 1822.888486, 1.00784 * 1822.888486,
		15.99491 * 1822.888486, 1.00784 * 1822.888486, 1.00784 * 1822.888486,
                15.99491 * 1822.888486, 1.00784 * 1822.888486, 1.00784 * 1822.888486,
                15.99491 * 1822.888486, 1.00784 * 1822.888486, 1.00784 * 1822.888486,
                15.99491 * 1822.888486, 1.00784 * 1822.888486, 1.00784 * 1822.888486,
                15.99491 * 1822.888486, 1.00784 * 1822.888486, 1.00784 * 1822.888486,
                15.99491 * 1822.888486, 1.00784 * 1822.888486, 1.00784 * 1822.888486,
                15.99491 * 1822.888486, 1.00784 * 1822.888486, 1.00784 * 1822.888486,
                15.99491 * 1822.888486, 1.00784 * 1822.888486, 1.00784 * 1822.888486,
                15.99491 * 1822.888486, 1.00784 * 1822.888486, 1.00784 * 1822.888486,
                15.99491 * 1822.888486, 1.00784 * 1822.888486, 1.00784 * 1822.888486],
    :target_state_index => 0, # Still needed for parsing logic if excited states were ever used
    :charge => 0,
    :multiplicity => 1,

    # --- Settings for Direct ORCA Execution ---
    :orca_executable => "/opt/orca-6.0.0/orca", # Full path to the ORCA binary
    :explicit_pal1 => true, # Add "! PAL1" to input if not found? Recommended.

    # --- Step Directory Settings ---
    :step_dir_prefix => "orca_calc_step", # Prefix for subdirs in :output_dir
    :cleanup_step_files => false,       # Delete extra files (.out, .tmp etc) from step dir? Keep .inp, .engrad, .gbw

    # Output Settings
    :output_prefix => "md_run", # Prefix for final MD results (trajectory, energy plot)
    :save_freq => 1,
    :print_freq => 1,
)

settings[:dt] = settings[:dt_fs] * 41.341374575751

# --- End User Defined Settings ---

# Basic validation
try
    # ... (keep existing file checks for pos, vel, template) ...
    # Simply check if the file exists at the path
    if !isfile(settings[:orca_executable])
        error("ORCA executable file not found at path: $(settings[:orca_executable])")
    end
    if !isdir(settings[:output_dir])
        # Base output dir will be created by run_md.jl if needed
        # No need to check here, but ensure run_md.jl does `mkpath(cfg[:output_dir])` early
    end
catch e
    println("Configuration Error: ", e)
    rethrow(e)
end

end # module SimulationConfig
