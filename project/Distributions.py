import math

import numpy as np
import scipy.stats
from sympy.parsing.sympy_parser import parse_expr

from sympy import lambdify

# ------------------------------Position distributions------------------------------#


class Coordinate_dist_3D:
    """Base class for 3D position distributions. Subclasses must implement the generate() method."""

    def __init__(self, box_size, N):
        # Accept either [Len_X, Units_X, Len_Y, Units_Y, Len_Z, Units_Z] or [Len_X, Len_Y, Len_Z]
        if len(box_size) == 6:
            self.len_X = float(box_size[0])
            self.len_Y = float(box_size[2])
            self.len_Z = float(box_size[4])
        elif len(box_size) == 3:
            self.len_X = float(box_size[0])
            self.len_Y = float(box_size[1])
            self.len_Z = float(box_size[2])
        else:
            raise ValueError("box_size must be a length-3 numeric tuple/list or length-6 [len,unit,len,unit,len,unit]")
        self.N = N

    def generate(self):
        """Generate and return N position coordinates."""
        raise NotImplementedError("Subclasses must implement generate()")


class Lattice_s(Coordinate_dist_3D):
    """Generate N positions arranged in a regular lattice within the box. Tries to fit to the box dimensions as best as possible."""

    def __init__(self, box_size, N):
        self.dist_name = "Lattice"
        super().__init__(box_size, N)
        self.lattice_info = self._find_best_lattice()

    def _find_best_lattice(self):
        """Find the best integer lattice factorization of N."""
        factors = []
        for nx in range(1, int(self.N ** (1 / 3)) + 2):
            if self.N % nx == 0:
                N2 = self.N // nx
                for ny in range(1, int(math.sqrt(N2)) + 2):
                    if N2 % ny == 0:
                        nz = N2 // ny
                        factors.append((nx, ny, nz))

        if len(factors) > 0:
            # Choose the factorization that gives the most cubic lattice (minimizes max spacing between particles in any dimension)
            best_factor = None
            best_spacing = float("inf")
            for nx, ny, nz in factors:
                ax = self.len_X / nx
                ay = self.len_Y / ny
                az = self.len_Z / nz
                max_spacing = max(ax, ay, az)
                if max_spacing < best_spacing:
                    best_spacing = max_spacing
                    best_factor = (nx, ny, nz)

            nx, ny, nz = best_factor
            ax = self.len_X / nx
            ay = self.len_Y / ny
            az = self.len_Z / nz
            return {
                "fits": True,
                "nx": nx,
                "ny": ny,
                "nz": nz,
                "lattice_spacing": (ax, ay, az),
            }
        return {"fits": False}

    def generate(self):
        """Generate N positions on a regular lattice."""
        if not self.lattice_info["fits"]:
            raise ValueError(f"Cannot arrange {self.N} particles in a cubic lattice")

        nx, ny, nz = (
            self.lattice_info["nx"],
            self.lattice_info["ny"],
            self.lattice_info["nz"],
        )
        ax, ay, az = self.lattice_info["lattice_spacing"]

        positions = []
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    x = (i + 0.5) * ax
                    y = (j + 0.5) * ay
                    z = (k + 0.5) * az
                    positions.append([x, y, z])

        return np.array(positions)


class Uniform_prob_density_s(Coordinate_dist_3D):
    """Generate N positions with uniform probability distribution in the box."""

    def __init__(self, box_size, N):
        self.dist_name = "Uniform pdf"
        super().__init__(box_size, N)

    def generate(self):
        """Generate N positions uniformly distributed in the box."""
        positions = np.zeros((self.N, 3))
        positions[:, 0] = np.random.uniform(0, self.len_X, self.N)
        positions[:, 1] = np.random.uniform(0, self.len_Y, self.N)
        positions[:, 2] = np.random.uniform(0, self.len_Z, self.N)
        return positions


class Gradient_prob_density_s(Coordinate_dist_3D):
    """Generate N positions distributed in the box according to a custom probability density function that can be defined by the user."""

    def __init__(self, box_size, N, density_func=None, func_expr_str=None):
        self.dist_name = "Custom gradient pdf"
        super().__init__(box_size, N)
        # density_func should be a callable that takes (x, y, z) and returns probability density
        self.density_func = density_func if density_func else self._default_density
        self.func_expr_str = func_expr_str
        self._grid_shape = self._choose_grid_shape()

    def _default_density(self, x, y, z):
        return 1.0

    def _choose_grid_shape(self):
        """Choose a sampling grid fine enough to represent the density without excessive cost."""
        approx_cells_per_axis = int(np.clip(np.ceil(self.N ** (1 / 3)) * 4, 12, 48))
        return (approx_cells_per_axis, approx_cells_per_axis, approx_cells_per_axis)

    def _build_density_function(self, func_expr):
        density_func = lambdify(("x", "y", "z"), func_expr, modules=["numpy"])

        def wrapped(x, y, z):
            values = density_func(x, y, z)
            values = np.asarray(values, dtype=float)
            if values.shape == ():
                return np.full(np.shape(x), float(values), dtype=float)
            return values

        return wrapped

    def _evaluate_density(self, density_func, x, y, z):
        density = density_func(x, y, z)
        density = np.asarray(density, dtype=float)
        if density.shape != np.shape(x):
            density = np.broadcast_to(density, np.shape(x)).astype(float)
        if not np.all(np.isfinite(density)):
            raise ValueError("Density function must return finite numeric values across the box")
        if np.any(density < 0):
            raise ValueError("Density function must be non-negative across the box")
        return density

    def _resolve_density_expr(self):
        """Get and validate a symbolic density expression from file or user input."""
        while True:
            try:
                if self.func_expr_str:
                    func_str = self.func_expr_str
                else:
                    func_str = input("Enter a custom probability density function P(x,y,z){0 <= x,y,z <= box limits}:\n")

                func_expr = parse_expr(func_str, evaluate=False)

                test_density = float(func_expr.subs({"x": self.len_X / 2, "y": self.len_Y / 2, "z": self.len_Z / 2}))
                if not np.isfinite(test_density):
                    raise ValueError("Density function must return a finite numeric value")
                return func_expr
            except Exception as e:
                if self.func_expr_str:
                    raise ValueError(f"Invalid gradient expression '{self.func_expr_str}': {e}") from e
                print(f"Invalid function: {e}. Please try again.")

    def generate(self):
        """Generate positions from a custom density using a vectorized cell sampler."""
        func_expr = self._resolve_density_expr()
        density_func = self._build_density_function(func_expr)

        nx, ny, nz = self._grid_shape
        dx = self.len_X / nx
        dy = self.len_Y / ny
        dz = self.len_Z / nz

        x_centers = (np.arange(nx) + 0.5) * dx
        y_centers = (np.arange(ny) + 0.5) * dy
        z_centers = (np.arange(nz) + 0.5) * dz
        grid_x, grid_y, grid_z = np.meshgrid(x_centers, y_centers, z_centers, indexing="ij")

        density = self._evaluate_density(density_func, grid_x, grid_y, grid_z)
        cell_weights = density.ravel()
        total_weight = cell_weights.sum()
        if total_weight <= 0:
            raise ValueError("Density function must be positive somewhere inside the box")

        probabilities = cell_weights / total_weight
        chosen_cells = np.random.choice(cell_weights.size, size=self.N, p=probabilities)

        ix, iy, iz = np.unravel_index(chosen_cells, (nx, ny, nz))
        positions = np.empty((self.N, 3), dtype=float)
        positions[:, 0] = (ix + np.random.rand(self.N)) * dx
        positions[:, 1] = (iy + np.random.rand(self.N)) * dy
        positions[:, 2] = (iz + np.random.rand(self.N)) * dz
        return positions


# ------------------------------Velocity distributions------------------------------#


class Vector_dist_3D:
    """Base class for velocity distributions. Subclasses must implement the generate() method."""

    def __init__(self, N, mass, temp):
        self.N = N
        self.mass = mass
        self.temp = temp
        self.k_B = 1.380649e-23  # Boltzmann constant in J/K

    def generate(self):
        """Generate and return N velocity vectors."""
        raise NotImplementedError("Subclasses must implement generate()")

    def _random_unit_vectors(self, N):
        """Generate N random unit vectors uniformly distributed on a sphere."""
        # Use random approach
        phi = np.random.uniform(0, 2 * np.pi, N)
        cos_theta = np.random.uniform(-1, 1, N)
        sin_theta = np.sqrt(1 - cos_theta**2)

        vx = sin_theta * np.cos(phi)
        vy = sin_theta * np.sin(phi)
        vz = cos_theta

        return np.column_stack((vx, vy, vz))


class Uniform_identical_RMS_v(Vector_dist_3D):
    """N velocities all with the same RMS, random directions."""

    def __init__(self, N, mass, temp):
        self.dist_name = "Identical RMS"
        super().__init__(N, mass, temp)

    def generate(self):
        """Generate N velocities all with the same RMS, random directions."""
        # RMS velocity from temperature: v_rms = sqrt(3 * k_B * T / m)
        v_rms = np.sqrt(3 * self.k_B * self.temp / self.mass)

        # Get random unit vectors and scale by v_rms
        unit_vectors = self._random_unit_vectors(self.N)
        velocities = unit_vectors * v_rms

        return velocities


class Uniform_range_RMS_v(Vector_dist_3D):
    """N velocities with RMS uniformly distributed between 0 and v_max, random directions."""

    def __init__(self, N, mass, temp):
        self.dist_name = "Uniform RMS pdf"
        super().__init__(N, mass, temp)

    def generate(self):
        """Generate N velocities with RMS uniformly distributed between 0 and v_max, random directions."""
        # v_max based on temperature
        v_max = np.sqrt(3 * self.k_B * self.temp / self.mass)

        # Generate random RMS values between 0 and v_max
        v_rms_values = np.random.uniform(0, v_max, self.N)

        # Get random unit vectors
        unit_vectors = self._random_unit_vectors(self.N)

        # Scale each unit vector by its corresponding RMS
        velocities = unit_vectors * v_rms_values[:, np.newaxis]

        return velocities


class Maxwell_Boltzmann_v(Vector_dist_3D):
    """N velocities drawn from Maxwell-Boltzmann distribution."""

    def __init__(self, N, mass, temp):
        self.dist_name = "Maxwell-Boltzmann"
        super().__init__(N, mass, temp)

    def generate(self):
        """Generate N velocities from Maxwell-Boltzmann distribution."""
        # Standard deviation for each component: sigma = sqrt(k_B * T / m)
        sigma = np.sqrt(self.k_B * self.temp / self.mass)

        # Generate N velocities with 3 components each from normal distribution
        v_rms_values = scipy.stats.maxwell.rvs(scale=sigma, size=self.N)

        # Get random unit vectors
        unit_vectors = self._random_unit_vectors(self.N)

        # Scale each unit vector by its corresponding RMS
        velocities = unit_vectors * v_rms_values[:, np.newaxis]

        return velocities


# ------------------------------Angular Momentum distributions------------------------------#


class AngularMomentum_dist_3D:
    """Base class for angular momentum distributions. Subclasses must implement the generate() method."""

    def __init__(self, N, moments_of_inertia, temp):
        """
        Initialize angular momentum distribution.

        Args:
            N: Number of particles
            moments_of_inertia: [I_x, I_y, I_z] moments of inertia in kg*m²
            temp: Temperature in K
        """
        self.N = N
        self.moments_of_inertia = self._normalise_moments_of_inertia(moments_of_inertia)
        self.temp = max(float(temp), 0.0)
        self.k_B = 1.380649e-23  # Boltzmann constant in J/K

    def _normalise_moments_of_inertia(self, moments_of_inertia):
        """
        Ensure moments of inertia are non-negative and in SI units (kg*m^2).

        Some initialization files contain moments in Da*Angstrom^2. Molecular moments
        in SI are tiny, so values much larger than SI scale are converted using:
        1 Da*Angstrom^2 = 1.66053906660e-47 kg*m^2.
        """
        moi = np.array(moments_of_inertia, dtype=float)
        moi = np.where(np.isfinite(moi) & (moi > 0.0), moi, 0.0)

        # Heuristic for non-SI values (e.g. O(1) entries from Da*Angstrom^2 input)
        if np.any(moi > 1e-35):
            da_a2_to_kg_m2 = 1.66053906660e-47
            moi = moi * da_a2_to_kg_m2

        return moi

    def generate(self):
        """Generate and return N angular momentum vectors."""
        raise NotImplementedError("Subclasses must implement generate()")


class Isotropic_thermal_L(AngularMomentum_dist_3D):
    """
    Angular momenta from isotropic thermal distribution.

    Uses equipartition theorem: <L_i²> = I_i * k_B * T
    Each component is drawn from normal distribution with σ_i = √(I_i * k_B * T)
    """

    def __init__(self, N, moments_of_inertia, temp):
        self.dist_name = "Isotropic thermal"
        super().__init__(N, moments_of_inertia, temp)

    def generate(self):
        """
        Generate N angular momentum vectors from thermal isotropic distribution.

        For non-rotating particles (all moments of inertia = 0), returns zero vectors.
        For rotating particles, generates from Gaussian distribution with variances
        determined by equipartition theorem.

        Returns:
            Array of shape (N, 3) with angular momentum components
        """
        # Check if particle has rotational degrees of freedom
        if np.all(self.moments_of_inertia == 0):
            # Point particle, no rotation
            return np.zeros((self.N, 3))

        # Standard deviations for each component: σ_i = √(I_i * k_B * T)
        sigma = np.sqrt(self.moments_of_inertia * self.k_B * self.temp)

        # Generate N angular momenta with 3 components from normal distribution
        L_x = np.random.normal(0, sigma[0], self.N)
        L_y = np.random.normal(0, sigma[1], self.N)
        L_z = np.random.normal(0, sigma[2], self.N)

        angular_momenta = np.column_stack((L_x, L_y, L_z))
        return angular_momenta


class Rayleigh_magnitude_L(AngularMomentum_dist_3D):
    """
    Angular momenta with magnitude from Rayleigh distribution (like Maxwell-Boltzmann for angular momentum).

    Uses equipartition theorem to set expected kinetic energy:
    <KE_rot> = (1/2) * <L²/I> = (3/2) * k_B * T for 3 rotational DOF
    """

    def __init__(self, N, moments_of_inertia, temp):
        self.dist_name = "Rayleigh magnitude"
        super().__init__(N, moments_of_inertia, temp)

    def _random_unit_vectors(self, N):
        """Generate N random unit vectors uniformly distributed on a sphere."""
        phi = np.random.uniform(0, 2 * np.pi, N)
        cos_theta = np.random.uniform(-1, 1, N)
        sin_theta = np.sqrt(1 - cos_theta**2)

        vx = sin_theta * np.cos(phi)
        vy = sin_theta * np.sin(phi)
        vz = cos_theta

        return np.column_stack((vx, vy, vz))

    def generate(self):
        """
        Generate N angular momentum vectors with Rayleigh-distributed magnitudes.

        The overall magnitude is drawn from Rayleigh distribution, ensuring
        proper thermal energy distribution. Direction is isotropic.

        Returns:
            Array of shape (N, 3) with angular momentum components
        """
        # Check if particle has rotational degrees of freedom
        if np.all(self.moments_of_inertia == 0):
            # Point particle, no rotation
            return np.zeros((self.N, 3))

        # For a 3D rotor with effective moment of inertia, use equipartition:
        # We need to find scale parameter for Rayleigh distribution
        # For 3 rotational DOF: < L²/(2*I_eff) > = (3/2)*k_B*T
        # Using average of moments: I_eff = mean(I_x, I_y, I_z)
        I_eff = np.mean(self.moments_of_inertia[self.moments_of_inertia > 0]) if np.any(self.moments_of_inertia > 0) else 0.0

        # Scale parameter for Rayleigh: σ such that <L²> = 3*k_B*T*I_eff
        # For Rayleigh: <L²> = 2*σ²
        sigma = np.sqrt(1.5 * self.k_B * self.temp * I_eff)

        # Generate magnitudes from Rayleigh distribution
        L_magnitudes = scipy.stats.rayleigh.rvs(scale=sigma, size=self.N)

        # Get random unit vectors
        unit_vectors = self._random_unit_vectors(self.N)

        # Scale each unit vector by its corresponding magnitude
        angular_momenta = unit_vectors * L_magnitudes[:, np.newaxis]

        return angular_momenta


def generate_positions(Current_Particles, Box_Params):
    N, pos_dist_type = Current_Particles[1], Current_Particles[7]
    gradient_expr = Current_Particles[9] if len(Current_Particles) > 9 else None
    if pos_dist_type == 1:
        return Lattice_s(Box_Params, N).generate()
    elif pos_dist_type == 2:
        return Uniform_prob_density_s(Box_Params, N).generate()
    elif pos_dist_type == 3:
        return Gradient_prob_density_s(Box_Params, N, func_expr_str=gradient_expr).generate()
    else:
        raise ValueError(f"Invalid position distribution type: {pos_dist_type}")


def generate_velocities(Current_Particles):
    N, mass, temp, vel_dist_type = (
        Current_Particles[1],
        Current_Particles[2],
        Current_Particles[6],
        Current_Particles[8],
    )
    if vel_dist_type == 1:
        return Uniform_identical_RMS_v(N, mass, temp).generate()
    elif vel_dist_type == 2:
        return Uniform_range_RMS_v(N, mass, temp).generate()
    elif vel_dist_type == 3:
        return Maxwell_Boltzmann_v(N, mass, temp).generate()
    else:
        raise ValueError(f"Invalid velocity distribution type: {vel_dist_type}")


def generate_angular_momenta(Current_Particles, ang_mom_dist_type=1):
    """
    Generate angular momentum vectors for particles based on their degrees of freedom.

    Args:
        Current_Particles: List containing particle parameters including N, moments_of_inertia, and temp
        ang_mom_dist_type: Distribution type (1=Isotropic thermal, 2=Rayleigh magnitude)

    Returns:
        Array of shape (N, 3) with angular momentum vectors
    """
    N = Current_Particles[1]
    inertia_data = Current_Particles[4]
    moments_of_inertia = inertia_data["moments_of_inertia"]
    temp = Current_Particles[6]

    if ang_mom_dist_type == 1:
        return Isotropic_thermal_L(N, moments_of_inertia, temp).generate()
    elif ang_mom_dist_type == 2:
        return Rayleigh_magnitude_L(N, moments_of_inertia, temp).generate()
    else:
        raise ValueError(f"Invalid angular momentum distribution type: {ang_mom_dist_type}")
