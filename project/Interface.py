import os

import Distributions as DIST
import Initialiser as INIT
import Particles as PART
import Simulation as SIMU

# -----------------build simulation from the initialisation-----------------
Box_Params, Time_Params, Sep_Particles = INIT.Initialise()

All_particles = []  # list of all particle objects in the simulation, will be updated as particles are initialised and created during the simulation
print(f"different types of particles: {len(Sep_Particles)}")
print(Sep_Particles)
for p in range(len(Sep_Particles)):
    print(Sep_Particles[p])
    type = Sep_Particles[p][0]
    N = Sep_Particles[p][1]
    mass = Sep_Particles[p][2]
    degrees = Sep_Particles[p][3]
    inertia_data = Sep_Particles[p][4]  # metadata for inertia and vibrational modes
    radius = Sep_Particles[p][5]
    temp = Sep_Particles[p][6]
    pos_dist_type = Sep_Particles[p][7]
    vel_dist_type = Sep_Particles[p][8]
    current_particle = [type, N, mass, degrees, inertia_data, radius, temp, pos_dist_type, vel_dist_type]
    positions = DIST.generate_positions(current_particle, Box_Params)
    velocities = DIST.generate_velocities(current_particle)
    for i in range(N):  # create particle objects and add them to the list of all particles
        particle_i = PART.particle(i, mass, degrees, radius, positions[i], velocities[i], type)
        All_particles.append(particle_i)


# -----------------run simulation-----------------
Len_X, Units_X, Len_Y, Units_Y, Len_Z, Units_Z = Box_Params
Time, Units_Time, dt, Units_dt, T_plt, Units_T_plt = Time_Params

# Use input file name as simulation name when file-based initialisation was chosen.
init_file = getattr(INIT, "last_initialisation_file", None)
if init_file:
    sim_name = os.path.basename(init_file)
    if sim_name.lower().endswith(".txt"):
        sim_name = sim_name[:-4]  # remove .txt extension
else:
    sim_name = f"Sim_{Len_X}{Units_X}_x_{Len_Y}{Units_Y}_x_{Len_Z}{Units_Z},_t={Time}{Units_Time},_dt={dt}{Units_dt}"

my_simulation = SIMU.simulation(sim_name, box_params=Box_Params, All_particles=All_particles, Time_Params=Time_Params)
print(f"Running simulation: {sim_name}")
my_simulation.run_simulation()
print("Simulation complete")
