import os
import time

import numpy as np

try:
    from . import Distributions as DIST
    from . import Initialiser as INIT
    from . import Particles as PART
    from . import Simulation as SIMU
except ImportError:
    import Distributions as DIST
    import Initialiser as INIT
    import Particles as PART
    import Simulation as SIMU

# -----------------build simulation from the initialisation-----------------
Box_Params, Time_Params, Sep_Particles, seed, init_file = INIT.Initialise()

# Set the random seed for reproducible distributions
np.random.seed(seed)

All_particles = []  # list of all particle objects in the simulation, will be updated as particles are initialised and created during the simulation
print(f"different types of particles: {len(Sep_Particles)}")
print(Sep_Particles)
for p in range(len(Sep_Particles)):
    print(Sep_Particles[p])
    particle_type = Sep_Particles[p][0]
    N = Sep_Particles[p][1]
    mass = Sep_Particles[p][2]
    degrees = Sep_Particles[p][3]
    inertia_data = Sep_Particles[p][4]  # metadata for inertia and vibrational modes
    radius = Sep_Particles[p][5]
    temp = Sep_Particles[p][6]
    pos_dist_type = Sep_Particles[p][7]
    vel_dist_type = Sep_Particles[p][8]
    gradient_expr = Sep_Particles[p][9] if len(Sep_Particles[p]) > 9 else None
    current_particle = [particle_type, N, mass, degrees, inertia_data, radius, temp, pos_dist_type, vel_dist_type, gradient_expr]
    positions = DIST.generate_positions(current_particle, Box_Params)
    velocities = DIST.generate_velocities(current_particle)

    # Generate angular momenta based on temperature, moments of inertia, and vibrational modes
    # Distribution type 1: Isotropic thermal (from equipartition theorem)
    angular_momenta = DIST.generate_angular_momenta(current_particle, ang_mom_dist_type=1)

    for i in range(N):  # create particle objects and add them to the list of all particles
        particle_i = PART.particle(i, mass, inertia_data, radius, positions[i], velocities[i], angular_momenta[i], particle_type)
        All_particles.append(particle_i)


# -----------------run simulation-----------------
Len_X, Units_X, Len_Y, Units_Y, Len_Z, Units_Z = Box_Params
Time, Units_Time, dt, Units_dt, T_plt, Units_T_plt = Time_Params

# Use input file name as simulation name when file-based initialisation was chosen.
print(f"Initialisation file: {init_file}")
# print(f"Initialisation file type: {type(init_file)}")
if init_file:
    sim_name = os.path.basename(init_file)
    if sim_name.lower().endswith(".txt"):
        sim_name = sim_name[:-4]  # remove .txt extension
else:
    sim_name = input("Enter a name for this simulation run: ").strip()
    if not sim_name:
        sim_name = "unnamed_simulation"

my_simulation = SIMU.simulation(sim_name, box_params=Box_Params, All_particles=All_particles, Time_Params=Time_Params)
print(f"Running simulation: {sim_name}")
start_time = time.perf_counter()
my_simulation.run_simulation()
elapsed_time = time.perf_counter() - start_time
print("Simulation complete")
print(f"Time taken: {elapsed_time:.3f} s")
