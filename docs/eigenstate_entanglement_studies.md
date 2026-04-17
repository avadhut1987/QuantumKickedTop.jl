# Eigenstate Entanglement

This module provides tools for computing **entanglement properties of Floquet eigenstates** in the kicked top system.

Eigenstate entanglement plays an important role in the study of **quantum chaos**, as chaotic systems typically produce highly entangled eigenstates.

The functionality in this module focuses on computing the **average bipartite entanglement entropy of Floquet eigenstates**.

---

# Floquet Eigenstates

The dynamics of the quantum kicked top are governed by a **Floquet operator**

$$
U
$$

which describes the unitary evolution of the system over one driving period.

The eigenvalue problem

$$
U |\psi_i\rangle = e^{i\phi_i} |\psi_i\rangle
$$

defines the **Floquet eigenstates** \( |\psi_i\rangle \).

These eigenstates form the basis for studying **stationary properties of the driven quantum system**.

---

# Bipartite Entanglement

To quantify entanglement, the system is partitioned into two subsystems:

- subsystem \(Q\)
- the remaining degrees of freedom

The reduced density matrix for subsystem \(Q\) is obtained by tracing out the rest of the system:

$$
\rho_Q = \mathrm{Tr}_{\bar Q} (|\psi\rangle \langle \psi|)
$$

where

- \( |\psi\rangle \) is a Floquet eigenstate
- \( \mathrm{Tr}_{\bar Q} \) denotes the partial trace over the complement of subsystem \(Q\).

---

# Von Neumann Entanglement Entropy

The entanglement entropy is measured using the **von Neumann entropy**

$$
S = -\mathrm{Tr}(\rho_Q \log \rho_Q)
$$

This quantity measures the degree of **quantum correlations between the two subsystems**.

For highly chaotic systems, eigenstates typically exhibit **large entanglement entropy**.

---

# Average Eigenstate Entropy

The main function implemented in this module computes the **average entropy across all Floquet eigenstates**.

Mathematically,

$$
\bar{S} =
\frac{1}{N+1}
\sum_{i=1}^{N+1}
S(\rho_Q^{(i)})
$$

where

- \(N+1\) is the Hilbert space dimension
- \(S(\rho_Q^{(i)})\) is the entropy of the \(i\)-th eigenstate.

The entropy is normalized by

$$
\log_2(Q+1)
$$

so that the value lies between **0 and 1**.

---

# Implemented Function

The module exports the function: average_entropy(N, U; Q)

Parameters:

- `N` : system size parameter
- `U` : Floquet operator
- `Q` : subsystem size (optional)

Default: Q = div(N - 1, 2)

which corresponds to a **balanced bipartition** of the system.

---

# Example Usage

```julia
using QuantumKickedTop
using QuantumKickedTop.EigenstateEntanglement

# system size
N = 100

# construct Floquet operator
U = floquet_operator(...)

# compute average eigenstate entropy
S_avg = average_entropy(N, U)

println(S_avg)
```

The returned value is the **normalized average eigenstate entanglement entropy**.

## Scientific Context

Eigenstate entanglement has been widely used as a diagnostic for **quantum chaos and thermalization**.

Chaotic quantum systems typically exhibit eigenstates that behave similarly to **random states**, leading to high entanglement entropy.

The kicked top provides a simple and well-studied model for investigating these phenomena.