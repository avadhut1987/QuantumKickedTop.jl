# Quantum Dynamics

This module implements the **quantum dynamics of the kicked top system**, including:

- construction of the **Floquet operator**
- generation of **quantum states (spin coherent states)**
- time evolution of quantum states

It forms the core of the quantum simulation framework in this package.

---

# Floquet Dynamics

The quantum kicked top is a **periodically driven system**, and its dynamics are governed by a **Floquet operator**.

## Floquet Operator

The evolution over one driving period is given by

$$
U = e^{-i\frac{k'}{2j}J_x^2} e^{-i\frac{k}{2j}J_z^2} e^{-ipJ_y}
$$

where

- $J_x, J_y, J_z$ are angular momentum operators  
- $k$ and $k'$ are the nonlinear kicking strengths  
- $p$ is the rotation parameter  
- $j$ is the spin quantum number  

The state evolves as

$$
|\psi_{n+1}\rangle = U |\psi_n\rangle
$$

---

## Double Kicked Top (DKT)

This package also supports the **Double Kicked Top (DKT)**, which extends the above dynamics by introducing **multiple kicks per period**.

- QKT → single kick  
- DKT → multiple kicks  

The standard Quantum Kicked Top is recovered as a **special case** of the DKT.

---

# FloquetSystem Module

The file

```
FloquetSystem.jl
```

constructs the Floquet operator for given parameters.

## What it does

- Builds angular momentum operators $J_x, J_y, J_z$
- Constructs exponential operators
- Combines them to form the unitary Floquet operator $U$

## Output

The result is a **unitary matrix** of dimension

$$
2j + 1
$$

which represents the quantum evolution operator.

---

# Quantum States (PhiStates Module)

The file

```
PhiStates.jl
```

implements the construction of **quantum states used as initial conditions**.

These states are crucial for:

- initializing quantum simulations  
- studying entanglement generation  
- connecting quantum dynamics with classical phase space  

---

## Spin Coherent States

The primary states used are **spin coherent states**

$$
|\theta, \phi\rangle
$$

which represent states localized around a point on the unit sphere.

They provide a direct connection between:

- classical phase space → $(\theta, \phi)$  
- quantum state → $ |\theta, \phi\rangle $

---

## Mathematical Construction

The coherent state is expanded in the angular momentum basis:

$$
|\theta, \phi\rangle =
\sum_{m=-j}^{j}
c_m(\theta, \phi)\, |j, m\rangle
$$

where

$$
c_m =
\binom{2j}{j+m}^{1/2}
\left(\cos\frac{\theta}{2}\right)^{j+m}
\left(\sin\frac{\theta}{2}\right)^{j-m}
e^{-i m \phi}
$$

---

## Implementation Details

The module performs the following steps:

1. Iterates over basis states $m = -j, \ldots, j$

2. Computes coefficients:

   - binomial factor  
   - trigonometric weights  
   - phase factor $e^{-i m \phi}$

3. Constructs the state vector:

```julia
ψ = [c_m for m in -j:j]
```

4. Ensures normalization

---

## Why Coherent States Matter

Spin coherent states:

- are **maximally localized** in phase space  
- mimic **classical initial conditions**  
- are ideal for studying **quantum–classical correspondence**  

Under chaotic dynamics, these states:

- spread over Hilbert space  
- generate entanglement  
- exhibit signatures of quantum chaos  

---

# Time Evolution

Given:

- initial state $ |\psi_0\rangle $
- Floquet operator $U$

the state after $n$ kicks is

$$
|\psi_n\rangle = U^n |\psi_0\rangle
$$

This enables the study of:

- dynamical evolution  
- entanglement growth  
- chaotic behavior  

---

# Example Usage

```julia
using QuantumKickedTop
using QuantumKickedTop.FloquetSystem
using QuantumKickedTop.PhiStates

# system parameters
j = 10
p = pi/2
k = 3.0
kp = 0.0

# initial condition (classical angles)
θ = 0.3
φ = 1.2

# construct Floquet operator
U = floquet(j, p, k, kp)

# construct a spin coherent state
b = SpinBasis(s)
psi0 = coherentspinstate(b, θ, φ) 

# evolve state n-times
psi = U^n * psi0
```

---

# Connection to Classical Dynamics

This module provides the bridge between:

- **classical phase space** → $(\theta, \phi)$  
- **quantum state** → coherent state  
- **quantum evolution** → Floquet dynamics  

Thus, it enables direct investigation of:

- quantum–classical correspondence  
- signatures of chaos  
- entanglement generation  

---

# Related Modules

This module interacts with:

- `classical_dynamics/` → classical trajectories  
- `eigenstate_entanglement_studies/` → eigenstate entropy  
- `quantum_information_measures/` → entropy and correlations  
- `visualization/` → plotting results  

---

# Scientific Context

The quantum kicked top is a fundamental model in:

- quantum chaos  
- Floquet systems  
- nonlinear quantum dynamics  

It provides a minimal system to study:

- transition from integrability to chaos  
- emergence of statistical behavior  
- entanglement generation in chaotic systems  

---

# References
- A. Seshadri, V. Madhok, and A. Lakshminarayan, *Tripartite
mutual information, entanglement, and scrambling in permutation symmetric systems with an application to quantum chaos*, Phys. Rev. E **98, 052205 (2018)**.

- A. V. Purohit and U. T. Bhosale, *Study of the double kicked top: A classical and quantum perspective*, Physical Review E **112, 014217 (2025)**.
