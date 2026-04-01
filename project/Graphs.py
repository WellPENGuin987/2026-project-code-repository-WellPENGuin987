# Graphs module for plotting various distributions at each epoch of the simulation
import os
import re

import matplotlib.pyplot as plt
import numpy as np


class graph:
    def __init__(self, title, x_axis, y_axis, time_step=None):
        self.title = title
        self.x_axis = x_axis
        self.y_axis = y_axis
        self.time_step = time_step

    def plot(self, data):
        print("Base graph class does not implement plot method")

    def _safe_filename(self, title):
        # Replace Windows-invalid filename characters and trim trailing spaces/dots.
        cleaned = re.sub(r'[<>:"/\\|?*]', "_", str(title))
        cleaned = cleaned.rstrip(" .")
        return cleaned if cleaned else "graph"

    def _save_or_show(self):
        plt.tight_layout()
        plt.grid(True, alpha=0.3)
        if hasattr(self, "save_dir") and self.save_dir:
            os.makedirs(self.save_dir, exist_ok=True)
            # Use only the descriptive part after " - " (sim_name is already encoded
            # in the directory path), keeping filenames short enough for Windows MAX_PATH.
            title_for_file = self.title.split(" - ", 1)[-1] if " - " in self.title else self.title
            safe_title = self._safe_filename(title_for_file)
            plt.savefig(os.path.join(self.save_dir, f"{safe_title}.png"), dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()

    def _x_label(self):
        return "Time (s)" if self.time_step else "Epoch"

    def _line_x_values(self, n_points):
        if self.time_step:
            return np.arange(n_points, dtype=float) * float(self.time_step)
        return np.arange(n_points, dtype=float)

    def _image_x_extent(self, n_points):
        # Set extent so image columns are centered on the sampled x values.
        if n_points <= 0:
            return [0.0, 1.0]
        if self.time_step:
            half_step = 0.5 * float(self.time_step)
            return [-half_step, (n_points - 1) * float(self.time_step) + half_step]
        return [-0.5, n_points - 0.5]


class position_distribution(graph):
    def __init__(self, title, x_axis, y_axis, time_step=None):
        super().__init__(title, x_axis, y_axis, time_step)

    def plot(self, data):
        if not data or len(data) < 2:
            return  # Need at least 2 epochs for time evolution
        # Compute radial distances from origin for each particle at each epoch
        distances = []
        epoch_means = []
        for epoch_positions in data:
            epoch_distances = [np.linalg.norm(pos) for pos in epoch_positions]
            distances.extend(epoch_distances)
            epoch_means.append(np.mean(epoch_distances))
        if not distances:
            return
        min_dist = min(distances)
        max_dist = max(distances)
        bins = np.linspace(min_dist, max_dist, 21)  # 20 bins
        # Compute histogram for each epoch
        histograms = []
        for epoch_positions in data:
            epoch_distances = [np.linalg.norm(pos) for pos in epoch_positions]
            hist, _ = np.histogram(epoch_distances, bins=bins)
            histograms.append(hist)
        histograms = np.array(histograms)
        # Plot as image: x=epoch, y=distance, color=frequency
        plt.figure(figsize=(10, 6))
        x_min, x_max = self._image_x_extent(len(data))
        im = plt.imshow(histograms.T, aspect="auto", origin="lower", extent=[x_min, x_max, min_dist, max_dist], cmap="viridis")
        cbar = plt.colorbar(im, label="Number of Particles")
        cbar.ax.tick_params(labelsize=10)
        plt.xlabel(self._x_label(), fontsize=12)
        plt.ylabel("Radial Distance from Origin (m)", fontsize=12)
        plt.title(self.title, fontsize=14, fontweight="bold")
        # Overlay mean line
        x_vals = self._line_x_values(len(data))
        # Extend mean line to the right edge of the graph
        x_vals_extended = list(x_vals) + [x_max]
        epoch_means_extended = epoch_means + [epoch_means[-1]]
        plt.plot(x_vals_extended, epoch_means_extended, "r-", linewidth=2, label="Mean Radial Distance")
        plt.legend(fontsize=10)
        plt.tick_params(axis="both", which="major", labelsize=10)
        self._save_or_show()


class velocity_distribution(graph):
    def __init__(self, title, x_axis, y_axis, time_step=None):
        super().__init__(title, x_axis, y_axis, time_step)

    def plot(self, data):
        if not data or len(data) < 2:
            return  # Need at least 2 epochs for time evolution
        # Compute speeds for each particle at each epoch
        speeds = []
        epoch_means = []
        for epoch_velocities in data:
            epoch_speeds = [np.linalg.norm(vel) for vel in epoch_velocities]
            speeds.extend(epoch_speeds)
            epoch_means.append(np.mean(epoch_speeds))
        if not speeds:
            return
        min_speed = min(speeds)
        max_speed = max(speeds)
        bins = np.linspace(min_speed, max_speed, 21)  # 20 bins
        # Compute histogram for each epoch
        histograms = []
        for epoch_velocities in data:
            epoch_speeds = [np.linalg.norm(vel) for vel in epoch_velocities]
            hist, _ = np.histogram(epoch_speeds, bins=bins)
            histograms.append(hist)
        histograms = np.array(histograms)
        # Plot as image: x=epoch, y=speed, color=frequency
        plt.figure(figsize=(10, 6))
        x_min, x_max = self._image_x_extent(len(data))
        im = plt.imshow(histograms.T, aspect="auto", origin="lower", extent=[x_min, x_max, min_speed, max_speed], cmap="viridis")
        cbar = plt.colorbar(im, label="Number of Particles")
        cbar.ax.tick_params(labelsize=10)
        plt.xlabel(self._x_label(), fontsize=12)
        plt.ylabel("Particle Speed (m/s)", fontsize=12)
        plt.title(self.title, fontsize=14, fontweight="bold")
        # Overlay mean line
        x_vals = self._line_x_values(len(data))
        # Extend mean line to the right edge of the graph
        x_vals_extended = list(x_vals) + [x_max]
        epoch_means_extended = epoch_means + [epoch_means[-1]]
        plt.plot(x_vals_extended, epoch_means_extended, "r-", linewidth=2, label="Mean Speed")
        plt.legend(fontsize=10)
        plt.tick_params(axis="both", which="major", labelsize=10)
        self._save_or_show()


class energy_distribution(graph):
    def __init__(self, title, x_axis, y_axis, time_step=None):
        super().__init__(title, x_axis, y_axis, time_step)

    def plot(self, data):
        if not data or len(data) < 2:
            return  # Need at least 2 epochs for time evolution
        # data is list of lists of energies
        all_energies = [e for epoch in data for e in epoch]
        if not all_energies:
            return
        min_energy = min(all_energies)
        max_energy = max(all_energies)
        bins = np.linspace(min_energy, max_energy, 21)  # 20 bins
        # Compute histogram for each epoch
        histograms = []
        for epoch_energies in data:
            hist, _ = np.histogram(epoch_energies, bins=bins)
            histograms.append(hist)
        histograms = np.array(histograms)
        # Plot as image: x=epoch, y=energy, color=frequency
        plt.figure(figsize=(10, 6))
        x_min, x_max = self._image_x_extent(len(data))
        im = plt.imshow(histograms.T, aspect="auto", origin="lower", extent=[x_min, x_max, min_energy, max_energy], cmap="viridis")
        cbar = plt.colorbar(im, label="Number of Particles")
        cbar.ax.tick_params(labelsize=10)
        plt.xlabel(self._x_label(), fontsize=12)
        plt.ylabel("Kinetic Energy per Particle (J)", fontsize=12)
        plt.title(self.title, fontsize=14, fontweight="bold")
        # Overlay mean line
        epoch_means = [np.mean(epoch_energies) for epoch_energies in data]
        x_vals = self._line_x_values(len(data))
        # Extend mean line to the right edge of the graph
        x_vals_extended = list(x_vals) + [x_max]
        epoch_means_extended = epoch_means + [epoch_means[-1]]
        plt.plot(x_vals_extended, epoch_means_extended, "r-", linewidth=2, label="Mean Energy")
        plt.legend(fontsize=10)
        plt.tick_params(axis="both", which="major", labelsize=10)
        self._save_or_show()


class temperature_distribution(graph):
    def __init__(self, title, x_axis, y_axis, time_step=None):
        super().__init__(title, x_axis, y_axis, time_step)

    def plot(self, data):
        if not data or len(data) < 2:
            return  # Need at least 2 epochs for time evolution
        # data is list of lists of temperatures
        all_temps = [t for epoch in data for t in epoch]
        if not all_temps:
            return
        min_temp = min(all_temps)
        max_temp = max(all_temps)
        bins = np.linspace(min_temp, max_temp, 21)  # 20 bins
        # Compute histogram for each epoch
        histograms = []
        for epoch_temps in data:
            hist, _ = np.histogram(epoch_temps, bins=bins)
            histograms.append(hist)
        histograms = np.array(histograms)
        # Plot as image: x=epoch, y=temperature, color=frequency
        plt.figure(figsize=(10, 6))
        x_min, x_max = self._image_x_extent(len(data))
        im = plt.imshow(histograms.T, aspect="auto", origin="lower", extent=[x_min, x_max, min_temp, max_temp], cmap="viridis")
        cbar = plt.colorbar(im, label="Number of Particles")
        cbar.ax.tick_params(labelsize=10)
        plt.xlabel(self._x_label(), fontsize=12)
        plt.ylabel("Temperature (K)", fontsize=12)
        plt.title(self.title, fontsize=14, fontweight="bold")
        # Overlay mean line
        epoch_means = [np.mean(epoch_temps) for epoch_temps in data]
        x_vals = self._line_x_values(len(data))
        # Extend mean line to the right edge of the graph
        x_vals_extended = list(x_vals) + [x_max]
        epoch_means_extended = epoch_means + [epoch_means[-1]]
        plt.plot(x_vals_extended, epoch_means_extended, "r-", linewidth=2, label="Mean Temperature")
        plt.legend(fontsize=10)
        plt.tick_params(axis="both", which="major", labelsize=10)
        self._save_or_show()


class angular_momentum_distribution(graph):
    def __init__(self, title, x_axis, y_axis, time_step=None):
        super().__init__(title, x_axis, y_axis, time_step)

        def plot(self, data):
            if not data:
                return  # No data to plot


class total_angular_momentum_squared_distribution(graph):
    def __init__(self, title, x_axis, y_axis, time_step=None):
        super().__init__(title, x_axis, y_axis, time_step)

    def plot(self, data):
        if not data:
            return

        # Total angular momentum at each recorded time point is the vector sum across particles.
        total_L_squared = []
        for epoch_L_vectors in data:
            if epoch_L_vectors is None or len(epoch_L_vectors) == 0:
                total_L_squared.append(0.0)
                continue
            L_total = np.sum(np.array(epoch_L_vectors, dtype=float), axis=0)
            total_L_squared.append(float(np.dot(L_total, L_total)))

        plt.figure(figsize=(10, 6))
        x_vals = self._line_x_values(len(total_L_squared))
        plt.plot(x_vals, total_L_squared, "m-", linewidth=2, marker="s", markersize=4, label=r"Total $|L_{tot}|^2$")
        plt.xlabel(self._x_label(), fontsize=12)
        plt.ylabel(r"Total Angular Momentum Squared $(\mathrm{kg}^2\,\mathrm{m}^4/\mathrm{s}^2)$", fontsize=12)
        plt.title(self.title, fontsize=14, fontweight="bold")
        plt.legend(fontsize=10)
        plt.tick_params(axis="both", which="major", labelsize=10)
        self._save_or_show()


class total_rotational_energy_distribution(graph):
    def __init__(self, title, x_axis, y_axis, time_step=None):
        super().__init__(title, x_axis, y_axis, time_step)

    def plot(self, data):
        if not data:
            return

        total_rotational_energy = [float(np.sum(epoch_rot)) for epoch_rot in data]
        plt.figure(figsize=(10, 6))
        x_vals = self._line_x_values(len(total_rotational_energy))
        plt.plot(x_vals, total_rotational_energy, "c-", linewidth=2.2, marker="D", markersize=4, label=r"Total Rotational Energy")
        plt.xlabel(self._x_label(), fontsize=12)
        plt.ylabel(r"Total Rotational Kinetic Energy (J)", fontsize=12)
        plt.title(self.title, fontsize=14, fontweight="bold")
        plt.legend(fontsize=10)
        plt.tick_params(axis="both", which="major", labelsize=10)
        self._save_or_show()


class total_vibrational_energy_distribution(graph):
    def __init__(self, title, x_axis, y_axis, time_step=None):
        super().__init__(title, x_axis, y_axis, time_step)

    def plot(self, data):
        if not data:
            return

        total_vibrational_energy = [float(np.sum(epoch_vib)) for epoch_vib in data]
        plt.figure(figsize=(10, 6))
        x_vals = self._line_x_values(len(total_vibrational_energy))
        plt.plot(x_vals, total_vibrational_energy, "y-", linewidth=2.2, marker="^", markersize=4, label=r"Total Vibrational Energy")
        plt.xlabel(self._x_label(), fontsize=12)
        plt.ylabel(r"Total Vibrational Energy (J)", fontsize=12)
        plt.title(self.title, fontsize=14, fontweight="bold")
        plt.legend(fontsize=10)
        plt.tick_params(axis="both", which="major", labelsize=10)
        self._save_or_show()


class total_kinetic_energy_components_distribution(graph):
    def __init__(self, title, x_axis, y_axis, time_step=None):
        super().__init__(title, x_axis, y_axis, time_step)

    def plot(self, translational_data, rotational_data, vibrational_data):
        if not translational_data or not rotational_data or not vibrational_data:
            return

        n_points = min(len(translational_data), len(rotational_data), len(vibrational_data))
        total_trans = [float(np.sum(translational_data[i])) for i in range(n_points)]
        total_rot = [float(np.sum(rotational_data[i])) for i in range(n_points)]
        total_vib = [float(np.sum(vibrational_data[i])) for i in range(n_points)]
        total_all = [total_trans[i] + total_rot[i] + total_vib[i] for i in range(n_points)]

        plt.figure(figsize=(10, 6))
        x_vals = self._line_x_values(n_points)
        plt.plot(x_vals, total_trans, "b-", linewidth=1.8, label=r"Translational ($E_{trans}$)", marker="o", markersize=3)
        plt.plot(x_vals, total_rot, "c-", linewidth=1.8, label=r"Rotational ($E_{rot}$)", marker="s", markersize=3)
        plt.plot(x_vals, total_vib, "y-", linewidth=1.8, label=r"Vibrational ($E_{vib}$)", marker="^", markersize=3)
        plt.plot(x_vals, total_all, "k-", linewidth=2.4, label=r"Total Energy", marker="D", markersize=3)
        plt.xlabel(self._x_label(), fontsize=12)
        plt.ylabel("Total System Energy (J)", fontsize=12)
        plt.title(self.title, fontsize=14, fontweight="bold")
        plt.legend(fontsize=10, loc="best")
        plt.tick_params(axis="both", which="major", labelsize=10)
        self._save_or_show()


class pressure_distribution(graph):
    def __init__(self, title, x_axis, y_axis, time_step=None):
        super().__init__(title, x_axis, y_axis, time_step)

    def plot(self, data):
        if not data:
            return  # No data to plot
        # Plot average pressures over time for ideal and impulse calculations
        avg_ideal = [np.mean([p[0] for p in epoch]) for epoch in data]
        avg_impulse = [np.mean([p[1] for p in epoch]) for epoch in data]
        plt.figure(figsize=(10, 6))
        x_vals = self._line_x_values(len(data))
        plt.plot(x_vals, avg_ideal, "r-", linewidth=2, label="Ideal Gas (PV=NkT)")
        plt.plot(x_vals, avg_impulse, "b-", linewidth=2, label="Impulse-Based (Wall Collisions)")
        plt.xlabel(self._x_label(), fontsize=12)
        plt.ylabel("System Pressure (Pa)", fontsize=12)
        plt.title(self.title, fontsize=14, fontweight="bold")
        plt.legend(fontsize=10, loc="best")
        plt.tick_params(axis="both", which="major", labelsize=10)
        self._save_or_show()


class collision_count(graph):
    def __init__(self, title, x_axis, y_axis, time_step=None):
        super().__init__(title, x_axis, y_axis, time_step)

    def plot(self, data):
        if not data:
            return
        pp_data = data["pp"]
        wall_data = data["wall"]
        if not pp_data:
            return
        pp_totals = [int(np.sum(c)) for c in pp_data]
        wall_totals = [int(np.sum(c)) for c in wall_data]
        total = [p + w for p, w in zip(pp_totals, wall_totals)]
        x_vals = self._line_x_values(len(pp_data))
        plt.figure(figsize=(10, 6))
        plt.plot(x_vals, total, "k-", linewidth=2, marker="o", markersize=4, label="Total")
        plt.plot(x_vals, pp_totals, "b-", linewidth=1.5, marker="s", markersize=4, label="Particle-Particle")
        plt.plot(x_vals, wall_totals, "r-", linewidth=1.5, marker="^", markersize=4, label="Particle-Wall")
        plt.xlabel(self._x_label(), fontsize=12)
        plt.ylabel("Collisions Per Sample Interval", fontsize=12)
        plt.title(self.title, fontsize=14, fontweight="bold")
        plt.legend(fontsize=10)
        plt.tick_params(axis="both", which="major", labelsize=10)
        self._save_or_show()


class cumulative_collision_count(graph):
    def __init__(self, title, x_axis, y_axis, time_step=None):
        super().__init__(title, x_axis, y_axis, time_step)

    def plot(self, data):
        if not data:
            return
        pp_data = data["pp"]
        wall_data = data["wall"]
        if not pp_data:
            return
        pp_totals = np.cumsum([int(np.sum(c)) for c in pp_data])
        wall_totals = np.cumsum([int(np.sum(c)) for c in wall_data])
        total = pp_totals + wall_totals
        x_vals = self._line_x_values(len(pp_data))
        plt.figure(figsize=(10, 6))
        plt.plot(x_vals, total, "k-", linewidth=2, marker="o", markersize=4, label="Total")
        plt.plot(x_vals, pp_totals, "b-", linewidth=1.5, marker="s", markersize=4, label="Particle-Particle")
        plt.plot(x_vals, wall_totals, "r-", linewidth=1.5, marker="^", markersize=4, label="Particle-Wall")
        plt.xlabel(self._x_label(), fontsize=12)
        plt.ylabel("Cumulative Collision Count", fontsize=12)
        plt.title(self.title, fontsize=14, fontweight="bold")
        plt.legend(fontsize=10)
        plt.tick_params(axis="both", which="major", labelsize=10)
        self._save_or_show()


class pairwise_distance(graph):
    def __init__(self, title, x_axis, y_axis, time_step=None):
        super().__init__(title, x_axis, y_axis, time_step)

    def plot(self, data):
        if not data or len(data) < 2:
            return  # Need at least 2 epochs for time evolution
        # data is a list of numpy arrays — compute min/max without concatenating
        # to avoid creating a large temporary array or a Python list of scalars.
        if not any(len(epoch) for epoch in data):
            return
        min_dist = float(min(np.min(epoch) for epoch in data if len(epoch)))
        max_dist = float(max(np.max(epoch) for epoch in data if len(epoch)))
        bins = np.linspace(min_dist, max_dist, 21)  # 20 bins
        # Compute histogram for each epoch
        histograms = []
        for epoch_distances in data:
            hist, _ = np.histogram(epoch_distances, bins=bins)
            histograms.append(hist)
        histograms = np.array(histograms)
        # Plot as image: x=epoch, y=distance, color=frequency
        plt.figure(figsize=(10, 6))
        x_min, x_max = self._image_x_extent(len(data))
        im = plt.imshow(histograms.T, aspect="auto", origin="lower", extent=[x_min, x_max, min_dist, max_dist], cmap="viridis")
        cbar = plt.colorbar(im, label="Number of Particle Pairs")
        cbar.ax.tick_params(labelsize=10)
        plt.xlabel(self._x_label(), fontsize=12)
        plt.ylabel("Pairwise Separation Distance (m)", fontsize=12)
        plt.title(self.title, fontsize=14, fontweight="bold")
        # Overlay mean line
        epoch_means = [np.mean(epoch_distances) for epoch_distances in data]
        x_vals = self._line_x_values(len(data))
        # Extend mean line to the right edge of the graph
        x_vals_extended = list(x_vals) + [x_max]
        epoch_means_extended = epoch_means + [epoch_means[-1]]
        plt.plot(x_vals_extended, epoch_means_extended, "r-", linewidth=2, label="Mean Pairwise Distance")
        plt.legend(fontsize=10)
        plt.tick_params(axis="both", which="major", labelsize=10)
        self._save_or_show()


class nearest_neighbor_distance(graph):
    def __init__(self, title, x_axis, y_axis, time_step=None):
        super().__init__(title, x_axis, y_axis, time_step)

    def plot(self, data):
        if not data or len(data) < 2:
            return
        # data is a list of 1-D arrays: nearest-neighbour distance per particle at each epoch
        all_distances = np.concatenate([np.asarray(d) for d in data])
        if all_distances.size == 0:
            return
        min_dist = float(all_distances.min())
        max_dist = float(all_distances.max())
        if min_dist == max_dist:
            return
        bins = np.linspace(min_dist, max_dist, 21)  # 20 bins
        histograms = []
        epoch_means = []
        for epoch_dists in data:
            epoch_arr = np.asarray(epoch_dists)
            hist, _ = np.histogram(epoch_arr, bins=bins)
            histograms.append(hist)
            epoch_means.append(float(np.mean(epoch_arr)))
        histograms = np.array(histograms)
        plt.figure(figsize=(10, 6))
        x_min, x_max = self._image_x_extent(len(data))
        im = plt.imshow(
            histograms.T,
            aspect="auto",
            origin="lower",
            extent=[x_min, x_max, min_dist, max_dist],
            cmap="viridis",
        )
        cbar = plt.colorbar(im, label="Number of Particles")
        cbar.ax.tick_params(labelsize=10)
        plt.xlabel(self._x_label(), fontsize=12)
        plt.ylabel("Nearest Neighbour Distance (m)", fontsize=12)
        plt.title(self.title, fontsize=14, fontweight="bold")
        x_vals = self._line_x_values(len(data))
        x_vals_extended = list(x_vals) + [x_max]
        epoch_means_extended = epoch_means + [epoch_means[-1]]
        plt.plot(x_vals_extended, epoch_means_extended, "r-", linewidth=2, label="Mean Nearest Neighbour Distance")
        plt.legend(fontsize=10)
        plt.tick_params(axis="both", which="major", labelsize=10)
        self._save_or_show()
