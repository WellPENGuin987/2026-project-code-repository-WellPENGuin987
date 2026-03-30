import Graphs


class simulation:
    def __init__(self, sim_name, box_params, All_particles, Time_Params):
        self.sim_name = sim_name
        self.box_params = box_params
        self.box_size = [box_params[0], box_params[2], box_params[4]]  # list of box dimensions in x, y, z directions
        self.All_particles = All_particles  # list of particle objects, will be updated during the simulation
        self.Time_Params = Time_Params  # list of time parameters, will be used to determine the total run time and time step of the simulation

    def plot_graphs(
        self, epoch, position_dist_data, RMS_vel_dist_data, energy_dist_data, temp_dist_data, pressure_dist_data, collision_count_data, pairwise_dist_data, com_dist_data
    ):
        # create graphs using the recorded data
        position_graph = Graphs.position_distribution(f"{self.sim_name}_position_distribution_epoch_{epoch}", x_axis="Particle Type", y_axis="Position")
        position_graph.plot(position_dist_data)
        vel_graph = Graphs.velocity_distribution(f"{self.sim_name}_velocity_distribution_epoch_{epoch}", x_axis="Particle Type", y_axis="RMS Velocity")
        vel_graph.plot(RMS_vel_dist_data)
        energy_graph = Graphs.energy_distribution(f"{self.sim_name}_energy_distribution_epoch_{epoch}", x_axis="Particle Type", y_axis="Kinetic Energy")
        energy_graph.plot(energy_dist_data)
        temp_graph = Graphs.temperature_distribution(f"{self.sim_name}_temperature_distribution_epoch_{epoch}", x_axis="Particle Type", y_axis="Temperature")
        temp_graph.plot(temp_dist_data)
        pressure_graph = Graphs.pressure_distribution(f"{self.sim_name}_pressure_distribution_epoch_{epoch}", x_axis="Particle Type", y_axis="Pressure")
        pressure_graph.plot(pressure_dist_data)
        collision_graph = Graphs.collision_count(f"{self.sim_name}_collision_count_epoch_{epoch}", x_axis="Time", y_axis="Number of Collisions")
        collision_graph.plot(collision_count_data)
        pairwise_graph = Graphs.pairwise_distance(f"{self.sim_name}_pairwise_distance_epoch_{epoch}", x_axis="Particle Type", y_axis="Pairwise Distance")
        pairwise_graph.plot(pairwise_dist_data)

    def run_simulation(self):
        Time, Units_Time, dt, Units_dt, T_plt, Units_T_plt = self.Time_Params
        # convert time parameters to seconds for easier calculations
        time_unit_conversion = {"ns": 1e-9, "mis": 1e-6, "ms": 1e-3, "s": 1}
        Time_s = Time * time_unit_conversion[Units_Time]
        dt_s = dt * time_unit_conversion[Units_dt]
        T_plt_s = T_plt * time_unit_conversion[Units_T_plt]
        epochs = int(Time_s / dt_s)  # number of time steps in the simulation

        # lists to store data for graphs
        position_dist_data = []  # list of positions of all particles at each time step for position distribution graph
        RMS_vel_dist_data = []  # list of RMS velocities of all particles at each time step for velocity distribution graph
        energy_dist_data = []  # list of kinetic energies of all particles at each time step for energy distribution graph
        temp_dist_data = []  # list of temperatures of all particles at each time step for
        pressure_dist_data = []  # list of pressures of all particles at each time step for pressure distribution graph
        collision_count_data = []  # list of number of collisions at each time step for collision count graph
        pairwise_dist_data = []  # list of pairwise distances between all particles at each time step for pairwise distance graph
        com_dist_data = []  # list of distances of all particles from the center of mass at each time step for center of mass distance graph

        # plot initial graphs at the start of the simulation
        self.plot_graphs(
            epoch=0,
            position_dist_data=position_dist_data,
            RMS_vel_dist_data=RMS_vel_dist_data,
            energy_dist_data=energy_dist_data,
            temp_dist_data=temp_dist_data,
            pressure_dist_data=pressure_dist_data,
            collision_count_data=collision_count_data,
            pairwise_dist_data=pairwise_dist_data,
            com_dist_data=com_dist_data,
        )

        for epoch in range(epochs):
            # update positions of all particles
            for particle in self.All_particles:
                particle.update_position(dt_s)
            # check for collisions and update velocities accordingly
            for particle in self.All_particles:
                particle.check_collisions(self.All_particles, self.box_size)

            for particle in self.All_particles:
                particle.update_velocity(dt_s)

            # record data for analysis
            if epoch % int(T_plt_s / dt_s) == 0:  # if it's time to record data for plots
                energy_dist_data.append([particle.kinetic_energy() for particle in self.All_particles])
                temp_dist_data.append([particle.temperature() for particle in self.All_particles])
                pressure_dist_data.append([particle.pressure() for particle in self.All_particles])
                collision_count_data.append([particle.check_collisions(self.All_particles, self.box_size) for particle in self.All_particles])
                pairwise_dist_data.append(
                    [
                        ((particle_i.position - particle_j.position) ** 2).sum() ** 0.5
                        for i, particle_i in enumerate(self.All_particles)
                        for j, particle_j in enumerate(self.All_particles)
                        if i < j
                    ]
                )
                com_dist_data.append([((particle.position - self.center_of_mass()) ** 2).sum() ** 0.5 for particle in self.All_particles])
                # plot graphs at this epoch
                self.plot_graphs(
                    epoch=epoch,
                    position_dist_data=position_dist_data,
                    RMS_vel_dist_data=RMS_vel_dist_data,
                    energy_dist_data=energy_dist_data,
                    temp_dist_data=temp_dist_data,
                    pressure_dist_data=pressure_dist_data,
                    collision_count_data=collision_count_data,
                    pairwise_dist_data=pairwise_dist_data,
                    com_dist_data=com_dist_data,
                )

        # plot final graphs at the end of the simulation
        # self.plot_graphs()
