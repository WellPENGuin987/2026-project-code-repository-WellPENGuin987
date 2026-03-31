import os

import Graphs
import numpy as np

try:
    from . import Particles as PART
except ImportError:
    import Particles as PART


class simulation:
    def __init__(self, sim_name, box_params, All_particles, Time_Params):
        self.sim_name = sim_name
        self.box_params = box_params
        # Accept either [Len_X, Units_X, Len_Y, Units_Y, Len_Z, Units_Z] or [Len_X, Len_Y, Len_Z]
        if len(box_params) == 6:
            self.box_size = [box_params[0], box_params[2], box_params[4]]  # list of numeric box dimensions in x, y, z directions
        elif len(box_params) == 3:
            self.box_size = [box_params[0], box_params[1], box_params[2]]
        else:
            raise ValueError("box_params must be length 3 or 6")
        self.All_particles = All_particles  # list of particle objects, will be updated during the simulation
        self.particle_array = PART.ParticleArray(All_particles)  # vectorized particle array for efficient numpy operations
        self.Time_Params = Time_Params  # list of time parameters, will be used to determine the total run time and time step of the simulation
        self.volume = self.box_size[0] * self.box_size[1] * self.box_size[2]  # volume of the simulation box
        self.total_area = 2 * (self.box_size[0] * self.box_size[1] + self.box_size[1] * self.box_size[2] + self.box_size[2] * self.box_size[0])  # total surface area of the box

    def center_of_mass(self):
        """
        Calculate the center of mass of all particles in the simulation.
        Returns the position vector of the center of mass.
        """
        total_mass = sum(particle.mass for particle in self.All_particles)
        if total_mass == 0:
            return np.array([0.0, 0.0, 0.0])
        com = sum(particle.mass * particle.position for particle in self.All_particles) / total_mass
        return com

    def plot_graphs(
        self,
        epoch,
        position_dist_data,
        RMS_vel_dist_data,
        energy_dist_data,
        temp_dist_data,
        pressure_dist_data,
        collision_count_data,
        pairwise_dist_data,
        angular_momentum_dist_data=None,
        rotational_energy_dist_data=None,
        vibrational_energy_dist_data=None,
    ):
        # create graphs using the recorded data
        save_dir = os.path.join(os.path.dirname(__file__), "plotted_graphs", self.sim_name)
        os.makedirs(save_dir, exist_ok=True)
        time_step = self.Time_Params[4] * {"ns": 1e-9, "mis": 1e-6, "ms": 1e-3, "s": 1}[self.Time_Params[5]]  # T_plt_s
        position_graph = Graphs.position_distribution(
            f"{self.sim_name} - Particle Radial Distance Distribution Over Time", x_axis="Particle Type", y_axis="Position", time_step=time_step
        )
        position_graph.save_dir = save_dir
        position_graph.plot(position_dist_data)
        vel_graph = Graphs.velocity_distribution(f"{self.sim_name} - Particle Speed Distribution Over Time", x_axis="Particle Type", y_axis="RMS Velocity", time_step=time_step)
        vel_graph.save_dir = save_dir
        vel_graph.plot(RMS_vel_dist_data)
        energy_graph = Graphs.energy_distribution(
            f"{self.sim_name} - Kinetic Energy per Particle Distribution Over Time", x_axis="Particle Type", y_axis="Kinetic Energy", time_step=time_step
        )
        energy_graph.save_dir = save_dir
        energy_graph.plot(energy_dist_data)
        temp_graph = Graphs.temperature_distribution(f"{self.sim_name} - Temperature Distribution Over Time", x_axis="Particle Type", y_axis="Temperature", time_step=time_step)
        temp_graph.save_dir = save_dir
        temp_graph.plot(temp_dist_data)
        pressure_graph = Graphs.pressure_distribution(f"{self.sim_name} - System Pressure (Ideal vs Impulse-Based)", x_axis="Particle Type", y_axis="Pressure", time_step=time_step)
        pressure_graph.save_dir = save_dir
        pressure_graph.plot(pressure_dist_data)
        collision_graph = Graphs.collision_count(f"{self.sim_name} - Cumulative Collision Count Over Time", x_axis="Time", y_axis="Number of Collisions", time_step=time_step)
        collision_graph.save_dir = save_dir
        collision_graph.plot(collision_count_data)
        pairwise_graph = Graphs.pairwise_distance(
            f"{self.sim_name} - Pairwise Separation Distance Distribution", x_axis="Particle Type", y_axis="Pairwise Distance", time_step=time_step
        )
        pairwise_graph.save_dir = save_dir
        pairwise_graph.plot(pairwise_dist_data)

        # Plot angular momentum data if available
        if angular_momentum_dist_data is not None:
            print("Angular momentum magnitudes recorded during simulation:")
            print(f"Number of recorded time points: {len(angular_momentum_dist_data)}")
            if len(angular_momentum_dist_data) > 0:
                # Average angular momentum magnitude over all particles at each time point
                avg_L_per_timestep = [np.mean([np.linalg.norm(p_L) for p_L in particles_L]) for particles_L in angular_momentum_dist_data]
                print(f"Mean |L| range: {np.min(avg_L_per_timestep):.6e} to {np.max(avg_L_per_timestep):.6e}")

                total_L2_graph = Graphs.total_angular_momentum_squared_distribution(
                    f"{self.sim_name} - Total System Angular Momentum Squared (Conservation Check)",
                    x_axis="Time",
                    y_axis="|L_tot|^2",
                    time_step=time_step,
                )
                total_L2_graph.save_dir = save_dir
                total_L2_graph.plot(angular_momentum_dist_data)

        # Plot rotational kinetic energy data if available
        if rotational_energy_dist_data is not None:
            print("Rotational kinetic energies recorded during simulation:")
            print(f"Number of recorded time points: {len(rotational_energy_dist_data)}")
            if len(rotational_energy_dist_data) > 0:
                avg_rot_KE_per_timestep = [np.mean(rot_KE_list) for rot_KE_list in rotational_energy_dist_data]
                print(f"Mean rot KE range: {np.min(avg_rot_KE_per_timestep):.6e} to {np.max(avg_rot_KE_per_timestep):.6e}")

                total_rot_graph = Graphs.total_rotational_energy_distribution(
                    f"{self.sim_name} - Total Rotational Kinetic Energy Over Time",
                    x_axis="Time",
                    y_axis="E_rot",
                    time_step=time_step,
                )
                total_rot_graph.save_dir = save_dir
                total_rot_graph.plot(rotational_energy_dist_data)

        if vibrational_energy_dist_data is not None and len(vibrational_energy_dist_data) > 0:
            max_vib_energy = max(float(np.sum(epoch_vib)) for epoch_vib in vibrational_energy_dist_data)
            if max_vib_energy > 0.0:
                total_vib_graph = Graphs.total_vibrational_energy_distribution(
                    f"{self.sim_name} - Total Vibrational Energy Over Time",
                    x_axis="Time",
                    y_axis="E_vib",
                    time_step=time_step,
                )
                total_vib_graph.save_dir = save_dir
                total_vib_graph.plot(vibrational_energy_dist_data)

        if (
            rotational_energy_dist_data is not None
            and vibrational_energy_dist_data is not None
            and len(rotational_energy_dist_data) > 0
            and len(vibrational_energy_dist_data) > 0
            and len(energy_dist_data) > 0
        ):
            total_components_graph = Graphs.total_kinetic_energy_components_distribution(
                f"{self.sim_name} - Total Energy Components: Translational, Rotational & Vibrational",
                x_axis="Time",
                y_axis="Energy",
                time_step=time_step,
            )
            total_components_graph.save_dir = save_dir
            total_components_graph.plot(energy_dist_data, rotational_energy_dist_data, vibrational_energy_dist_data)

    def _convert_time_params_to_seconds(self):
        Time, Units_Time, dt, Units_dt, T_plt, Units_T_plt = self.Time_Params
        time_unit_conversion = {"ns": 1e-9, "mis": 1e-6, "ms": 1e-3, "s": 1}
        Time_s = Time * time_unit_conversion[Units_Time]
        dt_s = dt * time_unit_conversion[Units_dt]
        T_plt_s = T_plt * time_unit_conversion[Units_T_plt]
        return Time_s, dt_s, T_plt_s

    def _initialize_data_buffers(self):
        return {
            "position_dist_data": [],
            "RMS_vel_dist_data": [],
            "energy_dist_data": [],
            "temp_dist_data": [],
            "pressure_dist_data": [],
            "collision_count_data": [],
            "pairwise_dist_data": [],
            "angular_momentum_dist_data": [],
            "rotational_energy_dist_data": [],
            "vibrational_energy_dist_data": [],
        }

    def _advance_particles_one_epoch(self, dt_s):
        # Use vectorized position update for efficiency
        self.particle_array.update_positions_vectorized(dt_s)
        self.particle_array.resolve_collisions_spatial(self.box_size)

        for particle in self.All_particles:
            particle.update_velocity(dt_s)

    def _compute_pairwise_distances(self):
        # Use vectorized computation for pairwise distances
        self.particle_array.sync_from_particles()
        return self.particle_array.compute_pairwise_distances_vectorized()

    def _record_epoch_data(self, data, dt_recording):
        # Sync particle array from individual particles for vectorized calculations
        self.particle_array.sync_from_particles()

        # Use vectorized calculations for efficiency
        kinetic_energies = self.particle_array.compute_kinetic_energies_vectorized()
        temperatures = self.particle_array.compute_temperatures_vectorized()
        ideal_pressures = self.particle_array.compute_ideal_pressures_vectorized(self.volume)
        impulse_pressures = self.particle_array.compute_impulse_pressures_vectorized(self.total_area, dt_recording)
        rotational_energies = self.particle_array.compute_rotational_kinetic_energies_vectorized()
        vibrational_energies = self.particle_array.compute_vibrational_energies_vectorized()

        data["energy_dist_data"].append(kinetic_energies.tolist())
        data["temp_dist_data"].append(temperatures.tolist())
        data["pressure_dist_data"].append([[ideal_pressures[i], impulse_pressures[i]] for i in range(len(self.All_particles))])
        data["collision_count_data"].append(self.particle_array.collision_counts.tolist())

        # Reset collision counts
        self.particle_array.reset_collision_counts_vectorized()

        data["pairwise_dist_data"].append(self._compute_pairwise_distances())
        data["angular_momentum_dist_data"].append([particle.angular_momentum.copy() for particle in self.All_particles])
        data["rotational_energy_dist_data"].append(rotational_energies.tolist())
        data["vibrational_energy_dist_data"].append(vibrational_energies.tolist())

    def _print_angular_momentum_summary(self, angular_momentum_dist_data, rotational_energy_dist_data):
        print("\n" + "=" * 60)
        print("ANGULAR MOMENTUM SUMMARY STATISTICS")
        print("=" * 60)
        if angular_momentum_dist_data:
            all_L_vectors = np.array([L_vec for timestep_L in angular_momentum_dist_data for L_vec in timestep_L])
            L_magnitudes = np.linalg.norm(all_L_vectors, axis=1)

            print("Angular momentum magnitude statistics:")
            print(f"  Mean:     {np.mean(L_magnitudes):.6e} kg·m²/s")
            print(f"  Std Dev:  {np.std(L_magnitudes):.6e} kg·m²/s")
            print(f"  Min:      {np.min(L_magnitudes):.6e} kg·m²/s")
            print(f"  Max:      {np.max(L_magnitudes):.6e} kg·m²/s")

            print("\nAngular momentum component statistics:")
            for i, axis in enumerate(["x", "y", "z"]):
                print(f"  L_{axis}: mean={np.mean(all_L_vectors[:, i]):.6e}, std={np.std(all_L_vectors[:, i]):.6e}")

        if rotational_energy_dist_data:
            all_rot_KE = np.concatenate(rotational_energy_dist_data)
            non_zero_rot_KE = all_rot_KE[all_rot_KE > 0]
            print("\nRotational kinetic energy statistics:")
            print(f"  Mean:     {np.mean(all_rot_KE):.6e} J")
            print(f"  Std Dev:  {np.std(all_rot_KE):.6e} J")
            if len(non_zero_rot_KE) > 0:
                print(f"  Min (>0): {np.min(non_zero_rot_KE):.6e} J")
                print(f"  Max:      {np.max(all_rot_KE):.6e} J")
        print("=" * 60 + "\n")

    def _flatten_epoch_data(self, dataset):
        if not dataset:
            return np.array([], dtype=float)
        return np.concatenate(dataset).astype(float)

    def _print_scalar_stats(self, title, values, unit):
        if values.size == 0:
            return
        print(title)
        print(f"  Mean:     {np.mean(values):.6e} {unit}")
        print(f"  Std Dev:  {np.std(values):.6e} {unit}")
        print(f"  Min:      {np.min(values):.6e} {unit}")
        print(f"  Max:      {np.max(values):.6e} {unit}")

    def _values_by_type(self, dataset, particle_types):
        if not dataset:
            return {}
        by_type = {}
        for epoch_values in dataset:
            for idx, value in enumerate(epoch_values):
                p_type = particle_types[idx]
                by_type.setdefault(p_type, []).append(float(value))
        return by_type

    def _print_per_type_summary(self, energy_dist_data, temp_dist_data, rotational_energy_dist_data, vibrational_energy_dist_data):
        if not (self.All_particles and energy_dist_data and temp_dist_data):
            return

        particle_types = [str(particle.type) for particle in self.All_particles]
        trans_by_type = self._values_by_type(energy_dist_data, particle_types)
        temp_by_type = self._values_by_type(temp_dist_data, particle_types)
        rot_by_type = self._values_by_type(rotational_energy_dist_data, particle_types) if rotational_energy_dist_data else {}
        vib_by_type = self._values_by_type(vibrational_energy_dist_data, particle_types) if vibrational_energy_dist_data else {}

        print("\nPer-particle-type summary:")
        for p_type in sorted(trans_by_type.keys()):
            trans_vals = np.array(trans_by_type.get(p_type, []), dtype=float)
            temp_vals = np.array(temp_by_type.get(p_type, []), dtype=float)
            rot_vals = np.array(rot_by_type.get(p_type, []), dtype=float)
            vib_vals = np.array(vib_by_type.get(p_type, []), dtype=float)

            print(f"\n  Type {p_type}:")
            if trans_vals.size > 0:
                print(f"    E_trans mean/std/min/max: {np.mean(trans_vals):.6e} / {np.std(trans_vals):.6e} / {np.min(trans_vals):.6e} / {np.max(trans_vals):.6e} J")
            if rot_vals.size > 0:
                print(f"    E_rot   mean/std/min/max: {np.mean(rot_vals):.6e} / {np.std(rot_vals):.6e} / {np.min(rot_vals):.6e} / {np.max(rot_vals):.6e} J")
            if vib_vals.size > 0 and np.max(vib_vals) > 0.0:
                print(f"    E_vib   mean/std/min/max: {np.mean(vib_vals):.6e} / {np.std(vib_vals):.6e} / {np.min(vib_vals):.6e} / {np.max(vib_vals):.6e} J")
            if trans_vals.size > 0 and rot_vals.size > 0 and vib_vals.size > 0:
                n_vals = min(len(trans_vals), len(rot_vals), len(vib_vals))
                total_vals = trans_vals[:n_vals] + rot_vals[:n_vals] + vib_vals[:n_vals]
                print(f"    E_total mean/std/min/max: {np.mean(total_vals):.6e} / {np.std(total_vals):.6e} / {np.min(total_vals):.6e} / {np.max(total_vals):.6e} J")
            if temp_vals.size > 0:
                print(f"    T       mean/std/min/max: {np.mean(temp_vals):.6e} / {np.std(temp_vals):.6e} / {np.min(temp_vals):.6e} / {np.max(temp_vals):.6e} K")

    def _print_energy_temperature_summary(self, energy_dist_data, temp_dist_data, rotational_energy_dist_data, vibrational_energy_dist_data):
        print("\n" + "=" * 60)
        print("ENERGY AND TEMPERATURE SUMMARY STATISTICS")
        print("=" * 60)

        all_trans_KE = self._flatten_epoch_data(energy_dist_data)
        all_rot_KE = self._flatten_epoch_data(rotational_energy_dist_data)
        all_vib_E = self._flatten_epoch_data(vibrational_energy_dist_data)
        all_temps = self._flatten_epoch_data(temp_dist_data)

        self._print_scalar_stats("Translational kinetic energy statistics:", all_trans_KE, "J")

        if all_rot_KE.size > 0:
            print()
            self._print_scalar_stats("Rotational kinetic energy statistics:", all_rot_KE, "J")

        if all_vib_E.size > 0 and np.max(all_vib_E) > 0.0:
            print()
            self._print_scalar_stats("Vibrational energy statistics:", all_vib_E, "J")

        if all_trans_KE.size > 0 and all_rot_KE.size > 0 and all_vib_E.size > 0:
            n_values = min(len(all_trans_KE), len(all_rot_KE), len(all_vib_E))
            total_E = all_trans_KE[:n_values] + all_rot_KE[:n_values] + all_vib_E[:n_values]
            print()
            self._print_scalar_stats("Total energy (trans + rot + vib) statistics:", total_E, "J")

        if all_temps.size > 0:
            print()
            self._print_scalar_stats("Temperature statistics:", all_temps, "K")

        self._print_per_type_summary(energy_dist_data, temp_dist_data, rotational_energy_dist_data, vibrational_energy_dist_data)

        print("=" * 60 + "\n")

    def run_simulation(self):
        Time_s, dt_s, T_plt_s = self._convert_time_params_to_seconds()
        dt_recording = T_plt_s
        epochs = int(Time_s / dt_s)
        record_interval = max(1, int(T_plt_s / dt_s))

        data = self._initialize_data_buffers()

        for epoch in range(epochs):
            print(f"Epoch {epoch + 1} / {epochs}")
            self._advance_particles_one_epoch(dt_s)

            if epoch % record_interval == 0:
                data["position_dist_data"].append([p.position.copy() for p in self.All_particles])
                data["RMS_vel_dist_data"].append([p.velocity.copy() for p in self.All_particles])
                self._record_epoch_data(data, dt_recording)

        self._print_angular_momentum_summary(
            data["angular_momentum_dist_data"],
            data["rotational_energy_dist_data"],
        )
        self._print_energy_temperature_summary(
            data["energy_dist_data"],
            data["temp_dist_data"],
            data["rotational_energy_dist_data"],
            data["vibrational_energy_dist_data"],
        )

        # Generate final graphs showing the complete time evolution
        self.plot_graphs(
            epoch=epochs,
            position_dist_data=data["position_dist_data"],
            RMS_vel_dist_data=data["RMS_vel_dist_data"],
            energy_dist_data=data["energy_dist_data"],
            temp_dist_data=data["temp_dist_data"],
            pressure_dist_data=data["pressure_dist_data"],
            collision_count_data=data["collision_count_data"],
            pairwise_dist_data=data["pairwise_dist_data"],
            angular_momentum_dist_data=data["angular_momentum_dist_data"],
            rotational_energy_dist_data=data["rotational_energy_dist_data"],
            vibrational_energy_dist_data=data["vibrational_energy_dist_data"],
        )
