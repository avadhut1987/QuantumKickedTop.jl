# Quantum Information Measures

This module implements **quantum information theoretic measures** used to analyze the dynamics of the kicked top system.

It provides tools to compute:

- reduced density matrices  
- entanglement entropy  
- linear entropy  
- quantum discord  
- concurrence  

These quantities are essential for studying **entanglement generation, correlations, and signatures of quantum chaos**.

---

# Density Matrix Formalism

A pure quantum state is represented as

$$
|\psi\rangle
$$

The corresponding density matrix is

$$
\rho = |\psi\rangle \langle \psi|
$$

To study entanglement, the system is partitioned into two subsystems:

- subsystem \(Q\)
- its complement \(\bar{Q}\)

---

# Reduced Density Matrix

The reduced density matrix is obtained by performing a **partial trace**:

$$
\rho_Q = \mathrm{Tr}_{\bar Q}(\rho)
$$

This operation removes degrees of freedom associated with the complementary subsystem.

---
## Implementation

The function `compute_rhoQ` computes the reduced density matrix by:

1. Reshaping the full state into a bipartite structure  
2. Performing a partial trace over the complementary subsystem  
3. Returning the reduced density matrix  

This is the **core operation** underlying all entanglement measures.

---

# Von Neumann Entropy

The entanglement entropy is defined as

$$
S = -\mathrm{Tr}(\rho_Q \log \rho_Q)
$$

It measures the degree of **quantum entanglement** between subsystems.

---

## Implementation

The function `compute_von_neumann_entropy` performs:

1. Eigenvalue decomposition of $\rho_Q$  
2. Computes  

$$
S = -\sum_i \lambda_i \log \lambda_i
$$

where $\lambda_i$ are the eigenvalues of $\rho_Q$.

---

# Linear Entropy

The linear entropy is a simpler measure of mixedness:

$$
S_L = 1 - \mathrm{Tr}(\rho_Q^2)
$$

---

## Implementation

The function `compute_linear_entropy` computes:

- $\mathrm{Tr}(\rho_Q^2)$  

and returns

$$
S_L = 1 - \mathrm{Tr}(\rho_Q^2)
$$

---

# Quantum Discord

Quantum discord captures **non-classical correlations beyond entanglement**.

It is defined as the difference between total and classical correlations.

---

## Implementation

The function `compute_quantum_discord` performs:

- construction of reduced density matrices  
- entropy evaluations  
- optimization over measurement bases  

This makes it **computationally more intensive** than entropy-based measures.

---

# Concurrence

Concurrence is a measure of **two-qubit entanglement**.

It is defined as

$$
C = \max(0, \lambda_1 - \lambda_2 - \lambda_3 - \lambda_4)
$$

---

## Implementation

The function `compute_concurrence` performs:

1. Construction of the spin-flipped density matrix  
2. Eigenvalue computation  
3. Evaluation of the concurrence formula  
---
# Module Structure

The functionality is implemented across the directory

`src/quantum_information_measures/`

including

- `QuantumUtils.jl`
- `QuantumUtilsDiscord.jl`

---

# Example Workflow

```julia
using QuantumKickedTop

# system parameters
j  = 10
p  = 1.0
k  = 3.0
kp = 0.5   # second kick (DKT parameter)

# construct Floquet operator (DKT form)
U = floquet(j, p, k, kp)

# compute eigenstates
eig = eigen(U.data)

ψ = eig.vectors[:, 1]

# subsystem parameters
N = 2*j
Q = div(N - 1, 2)

# reduced density matrix
ρQ = compute_rhoQ(ψ, N, Q)

# entropy measures
S  = compute_von_neumann_entropy(ρQ)
SL = compute_linear_entropy(ρQ)

println(S, SL)
```

---

# Computational Considerations

- Entropy calculations scale with **Hilbert space dimension \(2j+1\)**  
- Partial trace operations are **computationally expensive**  
- Quantum discord involves **optimization** and is significantly slower  

For large systems, it is recommended to use:

- multi-threading  
- parameter sampling  
- precomputed datasets  

---

# Scientific Context

These measures are central to the study of:

- quantum chaos  
- entanglement dynamics  
- thermalization in closed quantum systems  
- quantum correlations in driven systems  

In chaotic regimes, one typically observes:

- rapid growth of entropy  
- high eigenstate entanglement  
- near-random state behavior  

---

# References

- Nielsen & Chuang, *Quantum Computation and Quantum Information*

- A. V. Purohit and U. T. Bhosale,  
  *Study of the double kicked top: A classical and quantum perspective*, Physical Review E **112, 014217 (2025)**