# Angular Momentum Physics Reference

## Equipartition Theorem & Thermal Distribution

### Energy Distribution
For a system in thermal equilibrium at temperature T:
- Each quadratic degree of freedom contributes (1/2)k_B·T to mean energy
- 3 translational DOF: ⟨(1/2)m·v²⟩ = (3/2)k_B·T
- 3 rotational DOF: ⟨(1/2)I·ω²⟩ = (3/2)k_B·T (for spherically symmetric molecules)

### Angular Momentum from Temperature

**For a single rotational axis i:**
```
Rotational KE: KE_rot,i = (1/2)·L_i²/I_i

Equipartition: ⟨KE_rot,i⟩ = (1/2)k_B·T

Therefore: ⟨L_i²⟩ = I_i·k_B·T
```

**Standard deviation of each component:**
```
σ_i = √(I_i·k_B·T)
```

Each component L_x, L_y, L_z is independently drawn from:
```
L_i ~ N(0, σ_i²)
```

## Moment of Inertia

### For Rigid Bodies
```
I = ∫∫∫ r_⊥² dm
```

For principal axes (x, y, z):
```
I_x = moment of inertia about x-axis
I_y = moment of inertia about y-axis
I_z = moment of inertia about z-axis
```

### Example Values (in kg·m²)
- Diatomic molecules (linear): I ≈ 10⁻⁴⁶ to 10⁻⁴⁵ kg·m²
- Polyatomic molecules: I ≈ 10⁻⁴⁵ to 10⁻⁴⁴ kg·m²

## Degrees of Freedom

### Translation
- Always 3 DOF: motion in x, y, z directions

### Rotation
- **Linear molecules** (I_z ≈ 0): 2 DOF (rotation about 2 axes)
- **Non-linear molecules**: 3 DOF (rotation about 3 axes)

### Vibration
- Each vibrational mode adds 2 DOF (kinetic + potential energy)

### Total DOF Example
- **Point particle** (monoatomic): 3 (translation only)
- **Linear molecule** (diatomic): 5 (3 translation + 2 rotation)
- **Non-linear molecule** (H₂O): 6 (3 translation + 3 rotation)
- **With vibrations**: Add 2 per vibrational mode

## Implementation in Code

### Particle Degrees Configuration
```python
degrees = [3, inertia_data]  # inertia_data['degrees'] tells us total DOF
```

### Angular Momentum Magnitude
```python
|L| = √(L_x² + L_y² + L_z²)
```

### Rotational Kinetic Energy
```python
KE_rot = (1/2)·Σ(L_i²/I_i)  for i in {x,y,z}
```

### Total Kinetic Energy
```python
KE_total = KE_trans + KE_rot
         = (1/2)m·v² + (1/2)·Σ(L_i²/I_i)
```

### Temperature from Energy
```python
T = (2·KE_total)/(k_B·DOF)
```

## Collision Angular Momentum Transfer

### Elastic Collision
When two particles collide elastically:
1. **Momentum conservation**: m₁v₁ + m₂v₂ = m₁v₁' + m₂v₂'
2. **Energy conservation**: (1/2)m₁v₁² + (1/2)m₂v₂² = (1/2)m₁v₁'² + (1/2)m₂v₂'²

### Torque from Collision
The collision impulse J creates torque:
```
τ = r × J
```

Where r is the effective moment arm (distance from center to collision point).

For spherical particles at contact:
```
|r| ≈ radius
```

Angular momentum change:
```
ΔL = τ·Δt ≈ transfer_factor × |τ| × collision_impulse
```

## Numerical Notes

### Small Angle Approximation
For typical molecular collisions:
- Collision time: τ_c ≈ 10⁻¹⁵ s
- Angular momentum: L ≈ 10⁻⁴⁰ kg·m²/s
- Torque: τ ≈ 10⁻²⁵ N·m

### Common Temperatures
- **Room temperature**: T ≈ 300 K → k_B·T ≈ 4.1 × 10⁻²¹ J
- **High temperature**: T ≈ 1000 K → k_B·T ≈ 1.4 × 10⁻²⁰ J
- **Very high temp**: T ≈ 10000 K → k_B·T ≈ 1.4 × 10⁻¹⁹ J

### Resulting Angular Momenta
```
For I ≈ 10⁻⁴⁶ kg·m² and T ≈ 300 K:
σ = √(I·k_B·T) = √(10⁻⁴⁶ × 4.1 × 10⁻²¹) ≈ 6.4 × 10⁻³⁴ kg·m²/s
```

## References for Further Study

1. **Statistical Mechanics**
   - Equipartition theorem: Each quadratic term contributes (1/2)k_B·T
   - Maxwell-Boltzmann velocity distribution
   - Rayleigh-distributed speeds from Gaussian velocity components

2. **Rigid Body Dynamics**
   - Angular momentum: L = I·ω or L = r × p
   - Rotational kinetic energy: KE_rot = (1/2)·I·ω² = L²/(2I)
   - Torque: τ = dL/dt = r × F

3. **Collision Mechanics**
   - Coefficient of restitution
   - Impact parameter
   - Angular impulse and momentum transfer
