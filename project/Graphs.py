# Graphs module for plotting various distributions at each epoch of the simulation
import os

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

    def _save_or_show(self):
        if hasattr(self, "save_dir") and self.save_dir:
            os.makedirs(self.save_dir, exist_ok=True)
            plt.savefig(os.path.join(self.save_dir, f"{self.title}.png"))
            plt.close()
        else:
            plt.show()


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
        plt.figure()
        x_max = len(data) * self.time_step if self.time_step else len(data)
        plt.imshow(histograms.T, aspect="auto", origin="lower", extent=[0, x_max, min_dist, max_dist], cmap="viridis")
        plt.colorbar(label="Frequency")
        plt.xlabel("Time (s)" if self.time_step else "Epoch")
        plt.ylabel("Radial Distance")
        plt.title(self.title)
        # Overlay mean line
        x_vals = [i * self.time_step for i in range(len(data))] if self.time_step else range(len(data))
        plt.plot(x_vals, epoch_means, "r-", linewidth=2, label="Mean Radial Distance")
        plt.legend()
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
        plt.figure()
        x_max = len(data) * self.time_step if self.time_step else len(data)
        plt.imshow(histograms.T, aspect="auto", origin="lower", extent=[0, x_max, min_speed, max_speed], cmap="viridis")
        plt.colorbar(label="Frequency")
        plt.xlabel("Time (s)" if self.time_step else "Epoch")
        plt.ylabel("Speed")
        plt.title(self.title)
        # Overlay mean line
        x_vals = [i * self.time_step for i in range(len(data))] if self.time_step else range(len(data))
        plt.plot(x_vals, epoch_means, "r-", linewidth=2, label="Mean Speed")
        plt.legend()
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
        plt.figure()
        x_max = len(data) * self.time_step if self.time_step else len(data)
        plt.imshow(histograms.T, aspect="auto", origin="lower", extent=[0, x_max, min_energy, max_energy], cmap="viridis")
        plt.colorbar(label="Frequency")
        plt.xlabel("Time (s)" if self.time_step else "Epoch")
        plt.ylabel("Energy")
        plt.title(self.title)
        # Overlay mean line
        epoch_means = [np.mean(epoch_energies) for epoch_energies in data]
        x_vals = [i * self.time_step for i in range(len(data))] if self.time_step else range(len(data))
        plt.plot(x_vals, epoch_means, "r-", linewidth=2, label="Mean Energy")
        plt.legend()
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
        plt.figure()
        x_max = len(data) * self.time_step if self.time_step else len(data)
        plt.imshow(histograms.T, aspect="auto", origin="lower", extent=[0, x_max, min_temp, max_temp], cmap="viridis")
        plt.colorbar(label="Frequency")
        plt.xlabel("Time (s)" if self.time_step else "Epoch")
        plt.ylabel("Temperature")
        plt.title(self.title)
        # Overlay mean line
        epoch_means = [np.mean(epoch_temps) for epoch_temps in data]
        x_vals = [i * self.time_step for i in range(len(data))] if self.time_step else range(len(data))
        plt.plot(x_vals, epoch_means, "r-", linewidth=2, label="Mean Temperature")
        plt.legend()
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
        plt.figure()
        x_vals = [i * self.time_step for i in range(len(data))] if self.time_step else range(len(data))
        plt.plot(x_vals, avg_ideal, "r-", label="Ideal Gas Pressure")
        plt.plot(x_vals, avg_impulse, "b-", label="Impulse-Based Pressure")
        plt.xlabel("Time (s)" if self.time_step else "Epoch")
        plt.ylabel("Average Pressure")
        plt.title(self.title)
        plt.legend()
        self._save_or_show()


class collision_count(graph):
    def __init__(self, title, x_axis, y_axis, time_step=None):
        super().__init__(title, x_axis, y_axis, time_step)

    def plot(self, data):
        if not data:
            return  # No data to plot
        # Plot total collisions over time
        total_collisions = [sum(c) for c in data]
        plt.figure()
        x_vals = [i * self.time_step for i in range(len(data))] if self.time_step else range(len(data))
        plt.plot(x_vals, total_collisions)
        plt.xlabel("Time (s)" if self.time_step else "Time Step")
        plt.ylabel("Total Collisions")
        plt.title(self.title)
        self._save_or_show()


class pairwise_distance(graph):
    def __init__(self, title, x_axis, y_axis, time_step=None):
        super().__init__(title, x_axis, y_axis, time_step)

    def plot(self, data):
        if not data or len(data) < 2:
            return  # Need at least 2 epochs for time evolution
        # data is list of lists of pairwise distances
        all_distances = [d for epoch in data for d in epoch]
        if not all_distances:
            return
        min_dist = min(all_distances)
        max_dist = max(all_distances)
        bins = np.linspace(min_dist, max_dist, 21)  # 20 bins
        # Compute histogram for each epoch
        histograms = []
        for epoch_distances in data:
            hist, _ = np.histogram(epoch_distances, bins=bins)
            histograms.append(hist)
        histograms = np.array(histograms)
        # Plot as image: x=epoch, y=distance, color=frequency
        plt.figure()
        x_max = len(data) * self.time_step if self.time_step else len(data)
        plt.imshow(histograms.T, aspect="auto", origin="lower", extent=[0, x_max, min_dist, max_dist], cmap="viridis")
        plt.colorbar(label="Frequency")
        plt.xlabel("Time (s)" if self.time_step else "Epoch")
        plt.ylabel("Pairwise Distance")
        plt.title(self.title)
        # Overlay mean line
        epoch_means = [np.mean(epoch_distances) for epoch_distances in data]
        x_vals = [i * self.time_step for i in range(len(data))] if self.time_step else range(len(data))
        plt.plot(x_vals, epoch_means, "r-", linewidth=2, label="Mean Pairwise Distance")
        plt.legend()
        self._save_or_show()
