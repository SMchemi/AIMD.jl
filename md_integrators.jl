# md_integrators.jl
module MDIntegrators

export velocity_verlet_step!, calculate_kinetic_energy

"""
    velocity_verlet_step!(positions, velocities, forces, masses, dt)

Performs one step of the Velocity Verlet algorithm.
Updates positions and velocities in-place.
Assumes consistent atomic units.
positions, velocities, forces are [3 x N_atoms].
masses is [N_atoms].
"""
function velocity_verlet_step!(positions::Matrix{Float64},
                               velocities::Matrix{Float64},
                               forces::Matrix{Float64},
                               masses::Vector{Float64},
                               dt::Float64)
    # v(t + dt/2) = v(t) + 0.5 * a(t) * dt
    # Need acceleration a = F/m. Reshape masses for broadcasting over coordinates.
    mass_reshaped = reshape(masses, 1, length(masses))
    acceleration = forces ./ mass_reshaped # Element-wise division [3 x N] / [1 x N] -> [3 x N]
    velocities .+= 0.5 .* acceleration .* dt

    # r(t + dt) = r(t) + v(t + dt/2) * dt
    positions .+= velocities .* dt

    # NOTE: The second velocity update v(t + dt) = v(t + dt/2) + 0.5 * a(t + dt) * dt
    # requires the *new* forces a(t + dt), which must be calculated *after* the position update.
    # This step is completed in the main loop after the force calculation.
end

"""
    calculate_kinetic_energy(masses::Vector{Float64}, velocities::Matrix{Float64})

Calculates the total kinetic energy. Assumes atomic units.
Velocities [3 x N_atoms], masses [N_atoms]. Returns scalar energy in Hartree.
"""
function calculate_kinetic_energy(masses::Vector{Float64}, velocities::Matrix{Float64})
    # KE = 0.5 * sum(m_i * |v_i|^2)
    # |v_i|^2 = vx_i^2 + vy_i^2 + vz_i^2 = sum(velocities[:, i].^2)
    # Sum over coordinates first, then sum over atoms weighted by mass
    ke_per_atom = 0.5 .* masses .* vec(sum(velocities.^2, dims=1))
    return sum(ke_per_atom)
end

end # module MDIntegrators
