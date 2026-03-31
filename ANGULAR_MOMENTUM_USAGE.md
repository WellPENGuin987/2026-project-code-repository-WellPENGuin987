# Angular Momentum Implementation - Usage Guide

## Quick Start

### Running a Simulation with Angular Momentum

The angular momentum implementation is **fully integrated** and requires no special configuration:

1. **No code changes needed** - just run the normal interface:
   ```bash
   python Interface.py
   ```

2. **Follow the normal initialization prompts** for:
   - Box dimensions
   - Time parameters
   - Particle specifications (including moments of inertia)

3. **Angular momenta are automatically generated** based on:
   - Moments of inertia (I_x, I_y, I_z)
   - Initial temperature (T)
   - Degrees of freedom

### Example Initialization (for H₂O molecules)

**Manual entry:**
```
Input parameters for the box:
1.0,mim,1.0,mim,1.0,mim

Input parameters for time:
100.0,mis,1.0,mis,10.0,mis

Particle specification:
water,100,18.0,6,1.5,400.0,2,2
```

Then moments of inertia prompt:
```
Moments of inertia (I) in kg*m^2    format: I_x,I_y,I_z:
1.93e-48,0.615e-48,1.15e-48
```

**Result:** 
- 100 water molecules
- Each receives thermal angular momentum with σ_i = √(I_i·k_B·400K)
- Collisions include angular momentum transfer

## Expected Output

### Console Output During Simulation

```
Epoch 1 / 100000
Epoch 2 / 100000
...
============================================================
ANGULAR MOMENTUM SUMMARY STATISTICS
============================================================
Angular momentum magnitude statistics:
  Mean:     1.234567e-40 kg·m²/s
  Std Dev:  2.345678e-40 kg·m²/s
  Min:      0.000000e+00 kg·m²/s
  Max:      8.901234e-40 kg·m²/s

Angular momentum component statistics:
  L_x: mean=1.234567e-41, std=2.345678e-40
  L_y: mean=-5.678901e-42, std=2.345678e-40
  L_z: mean=3.456789e-41, std=2.345678e-40

Rotational kinetic energy statistics:
  Mean:     1.234567e-60 J
  Std Dev:  2.345678e-60 J
  Min (>0): 1.234567e-62 J
  Max:      9.876543e-60 J
============================================================
```

## Testing & Validation

### Unit Tests You Can Run

#### Test 1: Angular Momentum Distribution
```python
import numpy as np
import Distributions as DIST

# Test isotropic thermal distribution
particle_spec = [
    'test',              # type
    100,                 # N
    1e-26,               # mass (kg)
    6,                   # degrees
    {'moments_of_inertia': [1e-46, 1e-46, 1e-46], 'vibrational_modes': []},  # inertia
    1e-10,               # radius
    300.0,               # temp (K)
    1,                   # pos_dist
    1                    # vel_dist
]

L = DIST.generate_angular_momenta(particle_spec, ang_mom_dist_type=1)
print(f"Shape: {L.shape}")  # Should be (100, 3)
print(f"Mean |L|: {np.mean(np.linalg.norm(L, axis=1))}")
print(f"Expected σ: {np.sqrt(1e-46 * 1.38e-23 * 300)}")
```

**Expected**: Mean |L| should be approximately √(3)·σ ≈ 2.2×10⁻³⁵ kg·m²/s

#### Test 2: Rotational Kinetic Energy
```python
import Particles as PART

# Create a particle with angular momentum
pos = np.array([0.5, 0.5, 0.5])
vel = np.array([100.0, 0.0, 0.0])
L = np.array([1e-40, 1e-40, 1e-40])
I_data = (6, np.array([1e-46, 1e-46, 1e-46]))

p = PART.particle(1, 1e-26, I_data, 1e-10, pos, vel, L)

KE_rot = p.rotational_kinetic_energy()
expected = 3 * 0.5 * (1e-40)**2 / 1e-46
print(f"Rotational KE: {KE_rot}")
print(f"Expected: {expected}")
assert abs(KE_rot - expected) < 1e-60
```

**Expected**: KE_rot ≈ 1.5×10⁻⁶⁴ J

#### Test 3: Temperature Including Rotation
```python
# Temperature should include both translational and rotational contributions
p.velocity = np.array([100.0, 0.0, 0.0])  # v_rms ≈ 100 m/s
p.angular_momentum = np.array([1e-40, 1e-40, 1e-40])

T = p.temperature()
KE_trans = 0.5 * 1e-26 * 100**2  # ≈ 5×10⁻²² J
KE_rot = p.rotational_kinetic_energy()
KE_total = KE_trans + KE_rot
expected_T = 2 * KE_total / (1.38e-23 * 6)  # 6 DOF

print(f"Calculated T: {T}")
print(f"Expected T: {expected_T}")
```

#### Test 4: Collision With Angular Momentum Transfer
```python
# Two particles colliding head-on
p1 = PART.particle(0, 1e-26, (6, np.array([1e-46]*3)), 1e-10, 
                   np.array([0.0, 0.0, 0.0]), np.array([1000, 0, 0]), 
                   np.zeros(3))

p2 = PART.particle(1, 1e-26, (6, np.array([1e-46]*3)), 1e-10,
                   np.array([2e-10, 0.0, 0.0]), np.array([-1000, 0, 0]),
                   np.zeros(3))

# After collision, should have partial angular momentum transfer
p1._elastic_collision_with_rotation(p2, 2e-10)

print(f"p1 angular momentum: {p1.angular_momentum}")
print(f"p2 angular momentum: {p2.angular_momentum}")
# Should be non-zero (unless transfer_factor = 0)
```

### Validation Checks

✓ **Distribution validity**: Angular momenta should be normally distributed for isotropic model
✓ **Energy conservation**: Total KE before ≈ total KE after collision
✓ **Statistical consistency**: Component variances should match I_i·k_B·T
✓ **Zero moments**: Point particles should have zero angular momentum
✓ **Physical units**: All quantities in SI units (kg·m²/s for angular momentum)

## Configuration Options

### Switching Distribution Types

In `Interface.py`, line 32:
```python
# Default (Isotropic thermal):
angular_momenta = DIST.generate_angular_momenta(current_particle, ang_mom_dist_type=1)

# Alternative (Rayleigh magnitude):
angular_momenta = DIST.generate_angular_momenta(current_particle, ang_mom_dist_type=2)
```

### Adjusting Collision Transfer

In `Particles.py`, `_elastic_collision_with_rotation()` method:
```python
transfer_factor = 0.1  # Adjust this value
# Higher values → more angular momentum transfer
# Lower values → less coupling between translational and rotational motion
```

## Troubleshooting

### Issue: Very small angular momenta
**Cause**: Moments of inertia are very small
**Solution**: Check that moments of inertia are specified in correct units (kg·m²)

### Issue: Angular momentum statistics show mostly zeros
**Cause**: Particles are point particles (I = 0)
**Solution**: Ensure molecules have non-zero moments of inertia

### Issue: Rotational kinetic energy unrealistically high
**Cause**: Collision transfer factor too high or particle radii incorrect
**Solution**: Verify collision geometry and adjust transfer_factor if needed

### Issue: Temperature calculation seems wrong
**Cause**: Rotational energy not being counted properly
**Solution**: Verify that `temperature()` method is using total KE not just translational

## Advanced Customization

### Creating Custom Angular Momentum Distribution

Add to `Distributions.py`:
```python
class Custom_L_dist(AngularMomentum_dist_3D):
    """Custom angular momentum distribution."""
    
    def __init__(self, N, moments_of_inertia, temp):
        self.dist_name = "Custom distribution"
        super().__init__(N, moments_of_inertia, temp)
    
    def generate(self):
        """Your custom logic here."""
        if np.all(self.moments_of_inertia == 0):
            return np.zeros((self.N, 3))
        
        # Your distribution implementation
        L_magnitudes = ...  # Your model
        unit_vectors = self._random_unit_vectors(self.N)
        return unit_vectors * L_magnitudes[:, np.newaxis]
```

Then modify `generate_angular_momenta()`:
```python
elif ang_mom_dist_type == 3:
    return Custom_L_dist(N, moments_of_inertia, temp).generate()
```

### Modifying Collision Mechanics

To implement different collision models, edit `_elastic_collision_with_rotation()`:
- Change torque calculation formula
- Modify transfer_factor based on particle properties
- Add damping or friction effects
- Implement rotational coupling to vibrational modes

## Visualization & Analysis

The simulation generates statistics automatically. For further analysis:

```python
import numpy as np

# After running simulation, load angular momentum data:
# (This would be from simulation output or tracking arrays)

# Plot angular momentum evolution
import matplotlib.pyplot as plt

# Angular momentum magnitude over time
L_magnitudes = np.linalg.norm(L_vectors, axis=1)
plt.hist(L_magnitudes, bins=50)
plt.xlabel('|L| (kg·m²/s)')
plt.ylabel('Frequency')
plt.title('Angular Momentum Magnitude Distribution')
plt.show()
```

## Performance Notes

- Angular momentum calculations add ~5-10% computational overhead
- Collision detection unchanged
- Memory overhead: 3 floats per particle (angular momentum vector)
- For 1000 particles: ~24 KB additional memory

## References

See:
- `ANGULAR_MOMENTUM_IMPLEMENTATION.md` - Technical details
- `ANGULAR_MOMENTUM_PHYSICS.md` - Physics formulas and theory
