import math

import numpy as np
import scipy.stats
from sympy.parsing.sympy_parser import parse_expr

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

    def __init__(self, box_size, N, density_func=None):
        self.dist_name = "Custom gradient pdf"
        super().__init__(box_size, N)
        # density_func should be a callable that takes (x, y, z) and returns probability density
        self.density_func = density_func if density_func else self._default_density

    def generate(self):
        """Prompt user for custom density function."""

        while True:
            try:
                func_str = input("Enter a custom probability density function Ρ(x,y,z){0 ≤ x,y,z ≤ 1}:\n")
                # Parse the function string into a callable function
                func_expr = parse_expr(func_str, evaluate=False)

                # Test the function with a sample point to ensure it works and returns a numeric value
                test_density = np.float128(func_expr.subs({"x": self.len_X / 2, "y": self.len_Y / 2, "z": self.len_Z / 2}))
                if not isinstance(test_density, (int, float, np.float128)):
                    raise ValueError("Density function must return a numeric value")
                # normalise function over volume of box
                norm_const = 0
                print("Calculating normalization constant for the custom density function...")
                for _ in range(10000):  # Monte Carlo integration to find normalization constant
                    x = np.random.uniform(0, self.len_X)
                    y = np.random.uniform(0, self.len_Y)
                    z = np.random.uniform(0, self.len_Z)
                    norm_const += np.float128(func_expr.subs({"x": x, "y": y, "z": z}))
                norm_const *= (self.len_X * self.len_Y * self.len_Z) / 10000
                if norm_const == 0:
                    raise ValueError("Density function cannot be zero everywhere")

                break
            except Exception as e:
                print(f"Invalid function: {e}. Please try again.")
        positions = np.zeros((self.N, 3))
        for i in range(self.N):
            while True:
                x, y, z = np.random.rand(3)
                if np.random.rand() < np.float128(func_expr.subs({"x": self.len_X / 2, "y": self.len_Y / 2, "z": self.len_Z / 2})) / norm_const:
                    positions[i] = [x * self.len_X, y * self.len_Y, z * self.len_Z]
                    break
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


def generate_positions(Current_Particles, Box_Params):
    print(Current_Particles)
    N, pos_dist_type = Current_Particles[1], Current_Particles[7]
    if pos_dist_type == 1:
        return Lattice_s(Box_Params, N).generate()
    elif pos_dist_type == 2:
        return Uniform_prob_density_s(Box_Params, N).generate()
    elif pos_dist_type == 3:
        return Gradient_prob_density_s(Box_Params, N).generate()
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
