# MosaiQ-Lindhard: Berezinian (Superdeterminant) Architecture
**Target:** Electron-Phonon Supergraphs in Real-Axis Keldysh Dynamics  
**Status:** Phase 2 (Si Phonon Stabilization)  
**Author:** In-Gee Kim (Bona Sapiens, Inc.)  
**Docs map:** [`README.md`](README.md) · Phase~1 parent: [`simulator_architecture.md`](simulator_architecture.md)

## 1. Physical Motivation
In the strong-coupling regime of semiconductor physics (e.g., Silicon), the standard Random-Phase Approximation (RPA) catastrophically breaks down due to the topological overcounting of composite ring diagrams. This manifests macroscopically as the uniform hardening of phonon spectra and the breakdown of the Born-Oppenheimer approximation. 

To resolve this without introducing ad-hoc phenomenological parameters, we treat the coupled electron-phonon system as a **Graded Supergraph**. The statistical sign difference between fermions (electrons) and bosons (phonons) is strictly managed by substituting the conventional matrix determinant with the **Berezinian (Superdeterminant)**.

## 2. Mathematical Formulation
The coupled response is described by an asymmetric supermatrix $\mathcal{M}$:
$$
\mathcal{M} = \begin{pmatrix} A_{ee} & B_{ep} \\ C_{pe} & D_{pp} \end{pmatrix}
$$
*   **$A_{ee}$ (Fermion Block):** Electron-electron interactions ($N_e \times N_e$). Contains the fully relativistic 4-component Dirac-Gaussian LCAO matrix.
*   **$D_{pp}$ (Boson Block):** Phonon-phonon interactions ($N_p \times N_p$). Represents the phonon dynamical matrix.
*   **$B_{ep}, C_{pe}$ (Cross Blocks):** Electron-phonon couplings ($N_e \times N_p$ and $N_p \times N_e$).

The Berezinian is defined via the Schur complement:
$$
\operatorname{sdet}(\mathcal{M}) = \frac{\det(A_{ee} - B_{ep} D_{pp}^{-1} C_{pe})}{\det(D_{pp})}
$$

## 3. Numerical Strategy (Zero-Dependency C++20)
Direct computation of $D_{pp}^{-1}$ is numerically unstable and computationally wasteful. The MosaiQ engine enforces the following strict numerical policies:

### 3.1. The $DX = C$ Strategy (No Explicit Inversion)
Instead of inverting $D_{pp}$, we solve the linear system:
$$ D_{pp} X = C_{pe} $$
Then, the Schur complement is constructed via a single matrix multiplication:
$$ \tilde{A} = A_{ee} - B_{ep} X $$
Since $D_{pp}$ (Phonon Dynamical Matrix) is inherently symmetric/Hermitian, we utilize **LDLT (or Cholesky) decomposition** for optimal performance and stability.

### 3.2. Retarded Regularization ($i\eta$) for Causality
To preserve real-axis Keldysh causality and prevent singular inversions near resonance poles, a small imaginary broadening factor ($i\eta$) is diagonally loaded onto $D_{pp}$ before decomposition:
$$ D_{pp} \rightarrow D_{pp} + i\eta \mathbf{I} $$

### 3.3. Memory Architecture: `std::span`
To adhere to the Zero-Dependency philosophy, we avoid external tensor libraries (like Eigen) for core logic. The supermatrix is stored contiguously in memory, and blocks ($A_{ee}, D_{pp}, B_{ep}, C_{pe}$) are accessed safely without deep copies using C++20 `std::span`.

## 4. Implementation Blueprint (Header Interfaces)

Live tree (`simulator/src/`), aligned to Phase~1 layout:

*   `simulator/src/core/SuperMatrix.hpp`: Defines `ElectronPhononSuperMatrix<T>` with $N_e \neq N_p$ dimension support and `std::span` block views.
*   `simulator/src/core/LinalgUtils.hpp`: Implements standalone LU/LDLT decomposition, forward/backward substitution, and block multiplication (`DX=C`).
*   `simulator/src/engine/Berezinian.hpp`: Orchestrates the Schur complement computation (`computeSchurComplement`) applying the $DX=C$ strategy and $i\eta$ regularization.

---
*Note: This theoretical baseline directly supports the Keldysh dynamics transition (MosaiQ-Keldysh) and serves as the numerical proof for the forthcoming Springer Monograph.*
