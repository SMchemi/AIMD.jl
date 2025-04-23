# file_parsers.jl
module FileParsers

using DelimitedFiles
using Printf

export read_xyz, read_velocities, parse_comment_line

"""
Parses the second line (comment/title) of an XYZ/DAT file
expecting a format like:
'# ORCA AIMD Velocity Step X, t=Y fs, E_Pot=Z Hartree, Unit is ...'
Returns a dictionary with parsed values or defaults if parsing fails.
"""
function parse_comment_line(line::String)
    metadata = Dict{Symbol, Any}(:time_fs => nothing, :E_pot_Hartree => nothing)
    try
        # Regex to find t=VALUE fs and E_Pot=VALUE Hartree
        m_time = match(r"t=([\d\.\-eE]+)\s*fs", line)
        m_epot = match(r"E_Pot=([\d\.\-eE]+)\s*Hartree", line)

        if !isnothing(m_time)
            metadata[:time_fs] = parse(Float64, m_time.captures[1])
        end
        if !isnothing(m_epot)
            metadata[:E_pot_Hartree] = parse(Float64, m_epot.captures[1])
        end
    catch e
        @warn "Could not parse metadata from comment line: '$line'. Error: $e"
    end
    return metadata
end


"""
Reads atom symbols and coordinates from a standard XYZ file.
Also parses the comment line for optional metadata (time, energy).
Returns atom_symbols (Vector{String}) and coordinates (Matrix{Float64}, 3xN).
"""
function read_xyz(filename::String)
    lines = readlines(filename)
    if length(lines) < 3
        error("XYZ file '$filename' is too short.")
    end

    try
        n_atoms = parse(Int, lines[1])
        if length(lines) < n_atoms + 2
             error("XYZ file '$filename' has inconsistent number of atoms. Expected $n_atoms atoms based on line 1, but found fewer lines.")
        end

        comment_line = lines[2]
        metadata = parse_comment_line(comment_line)
        # You can use metadata if needed, e.g., print it:
        # println("Metadata from $filename: Time=$(metadata[:time_fs]) fs, Epot=$(metadata[:E_pot_Hartree]) Hartree")

        atom_symbols = Vector{String}(undef, n_atoms)
        coordinates = Matrix{Float64}(undef, 3, n_atoms)

        for i in 1:n_atoms
            parts = split(strip(lines[i+2]))
            if length(parts) < 4
                error("Error parsing line $(i+2) in '$filename'. Expected Atom X Y Z, got: $(lines[i+2])")
            end
            atom_symbols[i] = parts[1]
            coordinates[1, i] = parse(Float64, parts[2]) # X
            coordinates[2, i] = parse(Float64, parts[3]) # Y
            coordinates[3, i] = parse(Float64, parts[4]) # Z
        end

        return atom_symbols, coordinates

    catch e
        println("Error reading XYZ file '$filename': $e")
        rethrow(e)
    end
end

"""
Reads atom symbols and velocities from a DAT file assumed to be in XYZ format.
Verifies the number of atoms against the expected 'n_atoms_expected'.
Returns velocities (Matrix{Float64}, 3xN). Atom symbols are read but discarded.
"""
function read_velocities(filename::String, n_atoms_expected::Int)
    lines = readlines(filename)
    if length(lines) < 3
        error("Velocity file '$filename' is too short.")
    end

    try
        n_atoms_file = parse(Int, lines[1])
        if n_atoms_file != n_atoms_expected
            error("Mismatch in number of atoms between position file ($n_atoms_expected) and velocity file '$filename' ($n_atoms_file).")
        end
         if length(lines) < n_atoms_file + 2
             error("Velocity file '$filename' has inconsistent number of atoms. Expected $n_atoms_file atoms based on line 1, but found fewer lines.")
        end

        comment_line = lines[2]
        metadata = parse_comment_line(comment_line)
        # Optional: Use metadata if needed

        velocities = Matrix{Float64}(undef, 3, n_atoms_file)

        for i in 1:n_atoms_file
            parts = split(strip(lines[i+2]))
             if length(parts) < 4
                error("Error parsing line $(i+2) in '$filename'. Expected Atom Vx Vy Vz, got: $(lines[i+2])")
            end
            # We don't store the atom symbol from the velocity file, assume order matches XYZ
            velocities[1, i] = parse(Float64, parts[2]) # Vx
            velocities[2, i] = parse(Float64, parts[3]) # Vy
            velocities[3, i] = parse(Float64, parts[4]) # Vz
        end

        return velocities

    catch e
        println("Error reading velocity file '$filename': $e")
        rethrow(e)
    end
end

end # module FileParsers
