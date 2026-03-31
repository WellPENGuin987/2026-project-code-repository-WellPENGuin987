import numpy as np


class particle:
    """
    class representing a particle in the simulation, with attributes such as mass, position, velocity,
    and methods to update its position and velocity based on collisions with other particles and the walls of the box
    """

    def __init__(self, index, mass, degrees, radius, position, velocity, angular_momentum=None, type=None):
        self.index = index
        self.mass = mass
        if isinstance(degrees, dict):
            self.degrees = int(degrees.get("degrees", 3))
            moi = degrees.get("moments_of_inertia", [0.0, 0.0, 0.0])
            vib_modes = degrees.get("vibrational_modes", [])
        else:
            self.degrees = int(degrees)
            moi = [0.0, 0.0, 0.0]
            vib_modes = []
        self.moments_of_inertia = np.array(moi, dtype=float)
        self.vibrational_modes = np.array(vib_modes, dtype=float)
        self.radius = radius  # radius of collision, used to determine when particles collide
        self.position = np.array(position, dtype=float)
        self.velocity = np.array(velocity, dtype=float)
        self.angular_momentum = np.array(angular_momentum if angular_momentum is not None else [0.0, 0.0, 0.0], dtype=float)
        self.type = type if type is not None else index
        self.impulse_accum = 0  # accumulated impulses from wall collisions
        self.collision_count = 0  # number of collisions since last reset

    def update_position(self, dt):
        self.position += self.velocity * dt
        return self.position

    def check_collisions(self, All_particles, box_size):
        # check for collisions with other particles
        for particle in All_particles:
            if particle.index != self.index:  # don't check for collision with itself
                distance = ((self.position - particle.position) ** 2).sum() ** 0.5
                if distance < self.radius + particle.radius:  # collision detected
                    # update velocities of both particles according to elastic collision equations
                    m1 = self.mass
                    m2 = particle.mass
                    v1 = self.velocity
                    v2 = particle.velocity
                    self.velocity = v1 - (2 * m2 / (m1 + m2)) * ((v1 - v2).dot(self.position - particle.position) / distance**2) * (self.position - particle.position)
                    particle.velocity = v2 - (2 * m1 / (m1 + m2)) * ((v2 - v1).dot(particle.position - self.position) / distance**2) * (particle.position - self.position)
                    self.collision_count += 1
                    particle.collision_count += 1
        # check for collisions with walls of the box
        epoch_impulses = []  # list of impulses imparted on the particle by the walls during this epoch, used to calculate temperature changes
        for i in range(3):  # check x, y, z directions
            if self.position[i] < 0 or self.position[i] > box_size[i]:  # collision with wall detected
                old_v = self.velocity[i]
                self.velocity[i] = -self.velocity[i]  # reverse velocity in that direction
                epoch_impulses.append(2 * self.mass * abs(old_v))  # impulse imparted on the particle by the wall, used to calculate temperature changes
                self.collision_count += 1
        self.impulse_accum += sum(epoch_impulses)
        return epoch_impulses

    def update_velocity(self, dt):
        # For now, no additional velocity updates
        pass

    def kinetic_energy(self):
        return 0.5 * self.mass * np.sum(self.velocity**2)

    def rotational_kinetic_energy(self):
        positive_moi = self.moments_of_inertia > 0
        if not np.any(positive_moi):
            return 0.0
        return 0.5 * np.sum((self.angular_momentum[positive_moi] ** 2) / self.moments_of_inertia[positive_moi])

    def vibrational_energy(self):
        if self.vibrational_modes.size == 0:
            return 0.0
        # Classical equipartition estimate: each vibrational mode contributes k_B*T.
        k_B = 1.380649e-23
        return float(self.vibrational_modes.size * k_B * self.temperature())

    def temperature(self):
        # Simplified, assuming 3 degrees of freedom
        k_B = 1.380649e-23
        return (2 / 3) * self.kinetic_energy() / (k_B * self.degrees / 3)

    def pressure_ideal(self, volume):
        # Ideal gas pressure contribution per particle
        return (2 / 3) * self.kinetic_energy() / volume

    def pressure_impulse(self, area, dt):
        # Pressure from wall impulses
        p = self.impulse_accum / (6 * area * dt)
        self.impulse_accum = 0
        return p


class ParticleArray:
    """
    Vectorized particle management class for efficient numpy operations on all particles.
    Stores all particle data in contiguous arrays for fast vectorized operations.
    """

    def __init__(self, All_particles):
        """Initialize from a list of particle objects."""
        self.n_particles = len(All_particles)
        self.All_particles = All_particles  # Keep reference for compatibility

        # Vectorize basic properties
        self.indices = np.array([p.index for p in All_particles], dtype=int)
        self.masses = np.array([p.mass for p in All_particles], dtype=float)
        self.degrees = np.array([p.degrees for p in All_particles], dtype=int)
        self.radii = np.array([p.radius for p in All_particles], dtype=float)
        # Particle types can be strings (e.g. "water"), so keep them as generic labels.
        self.types = np.array([p.type for p in All_particles], dtype=object)

        # 2D arrays for per-particle scalar data
        self.moments_of_inertia = np.array([p.moments_of_inertia for p in All_particles], dtype=float)  # (n_particles, 3)
        # Store vibrational_modes as object array to handle variable lengths
        self.vibrational_modes = np.array([p.vibrational_modes for p in All_particles], dtype=object)

        # Position, velocity, angular momentum are stored in particle objects but synced here
        self.positions = np.array([p.position for p in All_particles], dtype=float)  # (n_particles, 3)
        self.velocities = np.array([p.velocity for p in All_particles], dtype=float)  # (n_particles, 3)
        self.angular_momenta = np.array([p.angular_momentum for p in All_particles], dtype=float)  # (n_particles, 3)

        self.impulse_accums = np.array([p.impulse_accum for p in All_particles], dtype=float)
        self.collision_counts = np.array([p.collision_count for p in All_particles], dtype=int)
        self._cell_size = max(1e-12, 2.0 * float(np.max(self.radii)) if self.n_particles else 1.0)

    def _cell_index(self, position, box_size):
        """Return integer cell coordinates for a position within the simulation box."""
        bounded = np.minimum(np.maximum(position, 0.0), np.array(box_size, dtype=float) - 1e-12)
        return tuple((bounded / self._cell_size).astype(int))

    def _build_spatial_cells(self, box_size):
        """Build a mapping from cell index to particle indices for broad-phase collision pruning."""
        cells = {}
        for i in range(self.n_particles):
            cell = self._cell_index(self.positions[i], box_size)
            if cell not in cells:
                cells[cell] = []
            cells[cell].append(i)
        return cells

    def _neighbor_cells(self, cell):
        """Generate all 27 neighboring cells (including the cell itself)."""
        cx, cy, cz = cell
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    yield (cx + dx, cy + dy, cz + dz)

    def _resolve_pair_collision(self, i, j):
        """Resolve an elastic collision between particles i and j if they overlap."""
        delta = self.positions[i] - self.positions[j]
        distance_sq = float(np.dot(delta, delta))
        min_dist = self.radii[i] + self.radii[j]

        if distance_sq <= 1e-24 or distance_sq >= (min_dist * min_dist):
            return

        m1 = self.masses[i]
        m2 = self.masses[j]
        v1 = self.velocities[i].copy()
        v2 = self.velocities[j].copy()

        factor_1 = (2.0 * m2 / (m1 + m2)) * (np.dot(v1 - v2, delta) / distance_sq)
        factor_2 = (2.0 * m1 / (m1 + m2)) * (np.dot(v2 - v1, -delta) / distance_sq)

        self.velocities[i] = v1 - factor_1 * delta
        self.velocities[j] = v2 - factor_2 * (-delta)
        self.collision_counts[i] += 1
        self.collision_counts[j] += 1

    def _resolve_neighbor_collisions(self, cells):
        """Resolve all candidate pair collisions from neighboring cell lists."""
        checked_pairs = set()
        for cell, indices in cells.items():
            for neighbor in self._neighbor_cells(cell):
                neighbor_indices = cells.get(neighbor)
                if neighbor_indices is None:
                    continue
                for i in indices:
                    for j in neighbor_indices:
                        if i >= j:
                            continue
                        pair = (i, j)
                        if pair in checked_pairs:
                            continue
                        checked_pairs.add(pair)
                        self._resolve_pair_collision(i, j)

    def _apply_wall_collisions(self, box_size):
        """Handle wall collisions vectorized across each Cartesian axis."""
        box = np.array(box_size, dtype=float)
        for axis in range(3):
            lower_hit = self.positions[:, axis] < 0.0
            upper_hit = self.positions[:, axis] > box[axis]
            hit = lower_hit | upper_hit
            if not np.any(hit):
                continue
            old_v = self.velocities[hit, axis].copy()
            self.velocities[hit, axis] *= -1.0
            self.impulse_accums[hit] += 2.0 * self.masses[hit] * np.abs(old_v)
            self.collision_counts[hit] += 1

    def resolve_collisions_spatial(self, box_size):
        """Resolve collisions using a cell-list broad phase plus exact pair checks."""
        if self.n_particles <= 1:
            return

        cells = self._build_spatial_cells(box_size)
        self._resolve_neighbor_collisions(cells)
        self._apply_wall_collisions(box_size)

        self.sync_to_particles()

    def sync_from_particles(self):
        """Sync arrays from particle objects (for non-vectorized operations)."""
        self.positions = np.array([p.position for p in self.All_particles], dtype=float)
        self.velocities = np.array([p.velocity for p in self.All_particles], dtype=float)
        self.angular_momenta = np.array([p.angular_momentum for p in self.All_particles], dtype=float)
        self.impulse_accums = np.array([p.impulse_accum for p in self.All_particles], dtype=float)
        self.collision_counts = np.array([p.collision_count for p in self.All_particles], dtype=int)

    def sync_to_particles(self):
        """Sync arrays back to particle objects."""
        for i, p in enumerate(self.All_particles):
            p.position = self.positions[i].copy()
            p.velocity = self.velocities[i].copy()
            p.angular_momentum = self.angular_momenta[i].copy()
            p.impulse_accum = self.impulse_accums[i]
            p.collision_count = self.collision_counts[i]

    def update_positions_vectorized(self, dt):
        """Vectorized position update: x += v * dt."""
        self.positions += self.velocities * dt
        # Sync back to particle objects
        for i, p in enumerate(self.All_particles):
            p.position[:] = self.positions[i]

    def compute_kinetic_energies_vectorized(self):
        """Vectorized kinetic energy calculation for all particles."""
        # KE = 0.5 * m * v^2 = 0.5 * m * sum(v_i^2)
        return 0.5 * self.masses * np.sum(self.velocities**2, axis=1)

    def compute_temperatures_vectorized(self):
        """Vectorized temperature calculation for all particles."""
        k_B = 1.380649e-23
        kinetic_energies = self.compute_kinetic_energies_vectorized()
        return (2 / 3) * kinetic_energies / (k_B * self.degrees / 3)

    def compute_rotational_kinetic_energies_vectorized(self):
        """Vectorized rotational kinetic energy calculation."""
        # Only include non-zero moments of inertia
        rotational_ke = np.zeros(self.n_particles)
        for i in range(self.n_particles):
            positive_moi = self.moments_of_inertia[i] > 0
            if np.any(positive_moi):
                rotational_ke[i] = 0.5 * np.sum((self.angular_momenta[i][positive_moi] ** 2) / self.moments_of_inertia[i][positive_moi])
        return rotational_ke

    def compute_vibrational_energies_vectorized(self):
        """Vectorized vibrational energy calculation."""
        k_B = 1.380649e-23
        temperatures = self.compute_temperatures_vectorized()
        n_vib_modes = np.array([len(vm) for vm in self.vibrational_modes], dtype=int)
        return n_vib_modes * k_B * temperatures

    def compute_pairwise_distances_vectorized(self):
        """Vectorized pairwise distance calculation using broadcasting."""
        # Compute all pairwise distances efficiently using broadcasting
        # Shape: (n_particles, n_particles) -> distances[i, j] = distance from particle i to j
        p_expanded_i = self.positions[:, np.newaxis, :]  # (n, 1, 3)
        p_expanded_j = self.positions[np.newaxis, :, :]  # (1, n, 3)
        differences = p_expanded_i - p_expanded_j  # (n, n, 3)
        distances_sq = np.sum(differences**2, axis=2)  # (n, n)
        distances = np.sqrt(distances_sq)

        # Extract upper triangle (i < j) to match original behavior
        upper_tri_indices = np.triu_indices(self.n_particles, k=1)
        return distances[upper_tri_indices]

    def compute_ideal_pressures_vectorized(self, volume):
        """Vectorized ideal gas pressure calculation."""
        kinetic_energies = self.compute_kinetic_energies_vectorized()
        return (2 / 3) * kinetic_energies / volume

    def compute_impulse_pressures_vectorized(self, area, dt):
        """Vectorized impulse-based pressure calculation."""
        pressures = self.impulse_accums / (6 * area * dt)
        self.impulse_accums[:] = 0  # Reset after computing
        for i, p in enumerate(self.All_particles):
            p.impulse_accum = 0
        return pressures

    def reset_collision_counts_vectorized(self):
        """Reset all collision counts."""
        self.collision_counts[:] = 0
        for p in self.All_particles:
            p.collision_count = 0
