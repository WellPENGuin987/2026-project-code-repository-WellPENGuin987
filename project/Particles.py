import numpy as np


class particle:
    """
    class representing a particle in the simulation, with attributes such as mass, position, velocity,
    and methods to update its position and velocity based on collisions with other particles and the walls of the box
    """

    def __init__(self, index, mass, degrees, radius, position, velocity, type=None):
        self.index = index
        self.mass = mass
        self.degrees = degrees  # degrees of freedom, used to determine specific heat capacity and thus temperature changes
        self.radius = radius  # radius of collision, used to determine when particles collide
        self.position = position
        self.velocity = velocity
        self.type = type if type is not None else index

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
                    pass
        # check for collisions with walls of the box
        epoch_impulses = []  # list of impulses imparted on the particle by the walls during this epoch, used to calculate temperature changes
        for i in range(3):  # check x, y, z directions
            if self.position[i] < 0 or self.position[i] > box_size[i]:  # collision with wall detected
                self.velocity[i] = -self.velocity[i]  # reverse velocity in that direction
                epoch_impulses.append(2 * self.mass * self.velocity[i])  # impulse imparted on the particle by the wall, used to calculate temperature changes
        return epoch_impulses

    def update_velocity(self, dt):
        # For now, no additional velocity updates
        pass

    def kinetic_energy(self):
        return 0.5 * self.mass * np.sum(self.velocity ** 2)

    def temperature(self):
        # Simplified, assuming 3 degrees of freedom
        k_B = 1.380649e-23
        return (2 / 3) * self.kinetic_energy() / (k_B * self.degrees / 3)

    def pressure(self):
        # Placeholder
        return 0
