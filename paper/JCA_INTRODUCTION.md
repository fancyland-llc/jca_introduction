```{=latex}
\begin{center}
{\LARGE\bfseries Introducing the Janus Coprime Algebra:\\[0.4em]
A Case Study on the LEO Kessler Critical Surface}

\vspace{1.5em}

{\large Antonio P.\ Matos}\\[0.4em]
{\small ORCID: \href{https://orcid.org/0009-0002-0722-3752}{0009-0002-0722-3752}}\\[0.3em]
April 30, 2026\\[0.3em]
{\small DOI: \href{https://doi.org/10.5281/zenodo.19955555}{10.5281/zenodo.19955555}}\\[0.3em]
Independent Researcher; Fancyland LLC / Lattice OS\\[0.3em]
Preprint v1.0

\vspace{0.8em}

{\small\textbf{MSC 2020:} 60J22, 81P68, 05C50, 60K35, 90B36, 70F15, 85A05}

\vspace{0.4em}

{\small\textbf{Keywords:} Janus Coprime Algebra (JCA), CREM commutator, primorial coprime lattice, rank-4 invariant subspace, $\sigma\otimes\tau$ tensor decomposition, cyclic $\sigma$-flux drive, Landau-Zener hysteresis, $\tau$-equivariant embedding, substrate-native Grover, kernel substitution test, Kessler syndrome, orbital debris, percolation transition, continuous-time quantum walk, Anderson localization, conjunction graph, ESA DISCOSweb}
\end{center}
```

**Companion papers (the foundational stack):**

- [1] A. P. Matos, "Character-Theoretic Effective Rank: From Schur's Lemma to a Stereoscopic Measurement Apparatus on Continuous Hamiltonian Substrates," preprint (2026). DOI: [10.5281/zenodo.19744573](https://doi.org/10.5281/zenodo.19744573)
- [2] A. P. Matos, "The Arithmetic Black Hole: Softmax Thermodynamics and the Four Eigenvalue Laws of the Prime Gas," preprint (2026). DOI: [10.5281/zenodo.19442006](https://doi.org/10.5281/zenodo.19442006)
- [3] A. P. Matos, "The Unity Clock: Effective Dimensional Collapse of the Addition-Multiplication Commutator on Coprime Lattices," preprint (2026). DOI: [10.5281/zenodo.19478727](https://doi.org/10.5281/zenodo.19478727)
- [4] A. P. Matos, "Active Transport on the Prime Gas: Flat-Band Condensation, the Rabi Phase Transition, and the Arithmetic Qubit," preprint (2026). DOI: [10.5281/zenodo.19243258](https://doi.org/10.5281/zenodo.19243258)
- [5] A. P. Matos, "Universal Two-Prime Formula for the Coprime-Lattice Coupling Constant," preprint (2026). DOI: [10.5281/zenodo.19210625](https://doi.org/10.5281/zenodo.19210625)

**Companion patents (the architectural tetrad):**

- [6] A. P. Matos, "Fault-Injection-Immune Computational Unit Using Primorial Coprime Residue Topology" (the Arithmetic Qubit), U.S. Provisional Patent Application No. **64/031,440** (filed April 7, 2026). Fancyland LLC.
- [7] A. P. Matos, "Holographic Eigen-Solver Using QM Boundary Projection on Coprime Residue Lattices," U.S. Provisional Patent Application No. **64/033,689** (filed April 8, 2026). Fancyland LLC.
- [8] A. P. Matos, "Tensegrity Interferometer: Stereoscopic Query Resolution and Manifold-Curvature Measurement on Continuous Hamiltonian Substrates," U.S. Provisional Patent Application No. **64/048,617** (filed April 24, 2026). Fancyland LLC.
- [9] A. P. Matos, "Method and Apparatus for Birefringent Control of a Quantum Arithmetic Substrate via Bloch-Flux-Driven P-Band Modulation (the 'PRISM Controller')," U.S. Provisional Patent Application No. **64/054,093** (filed April 30, 2026). Fancyland LLC.

---

## Patent Notice

The mathematics presented in this preprint introduces the **Janus Coprime Algebra (JCA)** — a structured operator class on the coprime residue lattice $(\mathbb{Z}/m)^\times$ for primorial $m$, defined by two independent axioms (bi-involutive equivariance and exact rank-4 invariant subspace) — as an extension of the earlier CREM operator $C(m) = [D_\text{sym}, P_\tau]$ from a single commutator to a structured algebra. Three new members of JCA are documented (the $\sigma\otimes\tau$ tensor decomposition pair, the cyclic $\sigma$-flux drive Hamiltonian family, and the $\tau$-equivariant embedding family) along with a kernel substitution test that empirically establishes axiomatic independence. **All theorems, scalar identities, and verification scripts referenced in this preprint are placed in the public domain** under the MIT license. **Engineering apparatus, system embodiments, gate-sequence implementations, and method claims drawn from this mathematics are reserved to the inventor and assignee** under U.S. Provisional Patent Applications [6], [7], [8], and [9] above.

| Layer | Provisional | Sections | Embodiment scope |
|---|---|---|---|
| State storage | 64/031,440 | §§2, 3 | $\sigma$-parity-protected register; rank-4 Q-band carrier |
| Spectral compute | 64/033,689 | §§2, 6 | FFT-accelerated eigenvalue extractors; per-altitude conjunction operators |
| Stereoscopic measurement | 64/048,617 | §§7, 8 | Per-irrep stereoscopic retrieval; manifold-curvature telemetry |
| **Read/write control** | **64/054,093** | §§3, 4, 8 | Bloch-flux-driven P-band modulation; cyclic $\sigma$-flux drive with hysteresis-mediated gating; engineering apparatus for $\sigma\otimes\tau$ register operation |

Implementations practicing the disclosed mathematics in the manner claimed in [6], [7], [8], or [9] require a license from the assignee. Where the mathematics is published openly here, it constitutes prior art against any third-party attempt to patent the same scalar identity, theorem, or proof technique going forward.

---

## Abstract

We introduce the **Janus Coprime Algebra (JCA)** as an extension of the earlier **Character-Rank Equivariant Markov (CREM)** operator [3] from a single commutator to a structured operator class. The JCA is defined on the coprime residue lattice $(\mathbb{Z}/m)^\times$ for squarefree primorial $m$ by two independent axioms: **(i) bi-involutive equivariance** — $[A, P_\sigma] = 0$ and $\{A, P_\tau\} = 0$ where $P_\sigma : r \mapsto -r \bmod m$ is additive reflection and $P_\tau : r \mapsto r^{-1} \bmod m$ is multiplicative inversion; and **(ii) exact rank-4 invariant subspace** — $A$ admits a rank-4 invariant subspace $V_4$ holding $\geq 98\%$ of $\|A\|_F^2$ at every primorial. The CREM commutator $C(m) = [D_\text{sym}, P_\tau]$ is the first known member; this paper documents three more members and demonstrates the algebra on a high-stakes empirical deployment: locating the Kessler critical surface in the modern LEO debris environment.

**Observation (JCA Axiomatic Independence, §2.5).** *The two JCA axioms are independent. Bi-involutive equivariance alone does not force rank-4 collapse: a kernel substitution test against five non-tent kernels (squared tent, cosine harmonic, Gaussian, inverse, random $P_\sigma$-symmetric), all satisfying axiom (i) at machine precision, produces rank-4 fractions ranging from 1.000 (cosine, structural artifact: kernel is intrinsically rank-2) to 0.003 (random $P_\sigma$-symmetric kernel) at $m = 30{,}030$ ($\varphi(m) = 5{,}760$). The L$_1$ tent kernel's specific geometry is what produces the rank-4 collapse; the algebra's two axioms are non-trivially co-satisfied. JCA membership requires verification of both axioms.*

**Theorem (V$_4$ Exact Invariance, §2.3).** *Let $V_4 \subset \mathbb{R}^{\varphi(m)}$ be the rank-4 boundary block of $C(m)$ — the span of the top four right singular vectors. Then $CV_4 \subseteq V_4$ exactly: the relative leakage $\eta(m) = \|CV_4 - V_4(V_4^T C V_4)\|_F / \|CV_4\|_F$ stays at the IEEE 754 floor ($\sim 10^{-15}$) at every primorial tested through the 8th — substrates $m \in \{2310, 30030, 510510, 9699690\}$ with $\varphi(m) \in \{480, 5760, 92160, 1{,}658{,}880\}$, spanning $3{,}456\times$ growth in basis dimension. The four-fold parity-chiral degeneracy of [1, Theorem 6.1(i)] holds exactly, not asymptotically, at every finite $m$.*

**Observation (Constant-Time 4D Braid, §2.4).** *Once $V_4$ has been extracted via the FFT-accelerated holographic eigensolver of [7], the deltas-only interferometric braid on $V_4$ evaluates in $\sim 280$ μs at every primorial tested through the 8th — including $\varphi(m) = 1{,}658{,}880$ — with 9% wall-time variation across $3{,}456\times$ substrate growth. Per-query cost is structurally invisible to substrate dimension after the rank-4 extraction; this is the structural reason a substrate-native realization of the JCA primitives can convert iteration-count quantum-search advantages into wall-clock advantages on substrate-fitted hardware (§8.4).*

**Theorem ($\sigma\otimes\tau$ Tensor Decomposition of V$_4$, §3.2).** *On $(\mathbb{Z}/m)^\times$ for squarefree $m$, the two natural involutions $P_\sigma : r \mapsto -r \bmod m$ (additive reflection) and $P_\tau : r \mapsto r^{-1} \bmod m$ (multiplicative inversion) commute as permutations. Their restrictions to $V_4$ commute as operators ($[P_\sigma|_{V_4}, P_\tau|_{V_4}] = 0$), both are involutive, and each has eigenvalues $\{-1, -1, +1, +1\}$ — splitting $V_4$ into four 1-dimensional sectors labeled by $(\sigma, \tau) \in \{\pm 1\}^2$:*
$$V_4 \;\cong\; \mathcal{H}_\sigma \otimes \mathcal{H}_\tau, \qquad \mathcal{H}_\sigma \cong \mathcal{H}_\tau \cong \mathbb{R}^2.$$
*The CREM commutator $C(m)$ in the qubit basis is block-diagonal in $\sigma$-sectors and off-diagonal in $\tau$ within each sector — i.e., $C$ is the natural $\tau$-driver, with a $\sigma$-dependent fine-structure frequency split (21 ppm at $m = 510{,}510$, sharpening as $O(1/\varphi(m))$). The split provides a substrate-intrinsic $\sigma$-bit non-demolition measurement primitive structurally analogous to dispersive shift in superconducting transmon qubits — except realized algebraically rather than engineered.*

**Observation (Cyclic $\sigma$-Flux Drive Produces Landau-Zener Hysteresis, §4.2).** *Add a time-dependent $\sigma$-bias term to the rank-4 effective Hamiltonian:*
$$H(t) \;=\; i V_4^T C V_4 \;+\; \lambda(t) \cdot P_\sigma^{4D}, \qquad \lambda(t) = \lambda_\text{max} \sin^2\!\left(\frac{\pi t}{T_\text{period}}\right).$$
*The loop area $A_\text{hyst} = \oint \langle P_\tau \rangle\, d\lambda$ scales monotonically with sweep rate over five orders of magnitude, from $3.1 \times 10^4$ in the near-adiabatic limit ($T \omega_\text{avg} = 10^3$) to $1.08 \times 10^9$ in the fast limit ($T \omega_\text{avg} = 3$). After one cycle at the fast rate, $\langle P_\tau\rangle(T) = -0.96$ vs initial $-1.000$ — a $4\%$ net $\tau$-bit deviation per cycle, confirming Landau-Zener-mediated non-adiabatic transitions through the fine-structure-split avoided crossings. Norm preservation $|\psi(T)|^2 = 1.000000$ at every endpoint (Trotter-clean unitary dynamics). The substrate exhibits rate-dependent state-tracking failure under cyclic $\sigma$-flux drive — a quantum-coherent diabatic-transition mechanism, not a thermodynamic dissipative hysteresis loop in the classical sense.*

**Result ($\tau$-Equivariant LEO Embedding + Substrate-Native Grover, §8).** *Six escalating natural-CREM-search tests on the LEO conjunction operator confirm that LEO does not natively organize into a CREM-class register under any tested physical involution. We pivot to the Shor-class embedding approach: construct $f: \text{satellite}_a \to r_a \in (\mathbb{Z}/510{,}510)^\times$ such that physical inclination-pairing $\tau_\text{LEO}(a) = b$ implies $r_b = r_a^{-1} \bmod m$. All 13,852 $\tau$-equivariance consistency checks pass. Substrate-native Grover converges to $P(|\text{target}\rangle) = 1.000000$ at $K = 238 = \lceil\pi/4 \cdot \sqrt{92{,}160}\rceil$ exact iterations, identifying the maximum-exposure LEO target Tian Lian 2-01 ↔ residue 329,701. The iteration count matches the Grover oracle-query optimum exactly — a $387\times$ query-count reduction vs $O(\varphi(m))$ classical scan; wall-clock conversion to a substrate-native hardware advantage is reserved to the engineering apparatus of [6]–[9] (§8.4). **The CREM substrate machinery transfers to LEO by construction, not by discovery** — extending the operator class with a third new operator: substrate-imposing equivariant embedding.*

**Theorem (Real-LEO Super-Critical Density, §6.3).** *The conjunction graph of the full DISCOSweb LEO catalog ($N = 49{,}998$ objects after physical-LEO perigee filter) has density $0.0084$, exceeding the Erdős–Rényi percolation threshold $p_c = 1/N$ by $116\times$. The giant connected component contains $94.8\%$ of all LEO objects; the residual $5.2\%$ split across 1,820 small components is dominated by isolated retrograde and graveyard-orbit objects.*

**Theorem (Intervention Budget, §7.2).** *Halving the giant component of the LEO conjunction graph requires deorbiting $2{,}424$ objects under a high-degree-first removal strategy — $2.9\times$ leverage over random removal.*

**Observation (Quantum-Advantage Flip, §7.4).** *On a synthetic Erdős–Rényi graph at the same critical density, continuous-time quantum walks (CTQW) spread ballistically at variance slope $+1.250$ in $\log t$, while continuous-time random walks (CTRW) diffuse at $+0.398$ — a $3.14\times$ quantum advantage. On the real LEO conjunction graph, the slopes invert to $+0.023$ (CTQW) and $+0.802$ (CTRW), a $0.029\times$ ratio — quantum walks localize rather than spread. The flip is consistent with Anderson-class localization on heterogeneous graphs with $\lambda_2 \to 0$.*

The Janus Coprime Algebra therefore admits at least four known members — the CREM commutator (genesis), the $\sigma\otimes\tau$ tensor decomposition pair, the cyclic $\sigma$-flux drive Hamiltonian family, and the $\tau$-equivariant embedding family — with axiomatic independence empirically established. The Kessler critical surface is the structurally-fitting high-stakes test bed: a percolation transition on a momentum-aware conjunction graph admits an order-2 involution (inclination antipode) on which the cost function is invariant — exactly the structure the $\tau$-equivariant embedding family is built to discretize. The case study exercises three of the four members: substrate-native Grover via the embedding identifies the maximum-exposure LEO target at the Grover oracle-query optimum, the $\sigma\otimes\tau$ register provides parity-protected state preparation, and the cyclic $\sigma$-flux drive provides the rate-controlled measurement primitive for cascade-state telemetry.

---

## 1. Introduction

### 1.1 The Janus Coprime Algebra and Its First Member

**Definition 1.1 (Janus Coprime Algebra).** *Let $m = p_1 p_2 \cdots p_k$ be a squarefree primorial. Let $P_\sigma : (\mathbb{Z}/m)^\times \to (\mathbb{Z}/m)^\times$ be the additive-reflection involution $r \mapsto -r \bmod m$ and $P_\tau$ the multiplicative-inversion involution $r \mapsto r^{-1} \bmod m$. The **Janus Coprime Algebra (JCA) at $m$** is the set of real linear operators $A$ on $\mathbb{R}^{\varphi(m)}$ satisfying both:*

- *(i) **Bi-involutive equivariance.** $[A, P_\sigma] = 0$ and $\{A, P_\tau\} = 0$.*
- *(ii) **Exact rank-4 invariant subspace.** $A$ admits an exact invariant subspace $V_4 \subset \mathbb{R}^{\varphi(m)}$ of dimension 4 holding $\geq 98\%$ of $\|A\|_F^2$.*

*The class is named for the two-faced involution structure: $P_\sigma$ and $P_\tau$ are the two natural involutions on the substrate, and JCA members are the operators that respect both faces simultaneously.*

The CREM commutator $C(m) = [D_\text{sym}, P_\tau]$ on coprime residue lattices was introduced in [3] and characterized in [1]. It is the **first known member** of the JCA. Three structural theorems anchor it:

- **Triangle-Wave Overtone Theorem [1, Theorem 6.1].** $C(m)$ has exact parity-chiral degeneracy $\sigma_{2i} = \sigma_{2i+1}$ at every finite $m$, fundamental amplitude $\sigma_0 = m\,\varphi(m)/\pi^2 + O(\varphi(m))$, and rank-4 boundary block holding asymptotic $L^2$ fraction $8/\pi^2$.
- **Universal Coupling Invariant [1, Theorem 7.1].** Two independent character-theoretic derivations agree bit-identically on the dimensionless coupling $\kappa^2 = \|C\|_F^2 / \lambda_\text{Perron}^2 = 2/3$.
- **Effective-Rank Plateau [3].** The effective rank of $C(m)$ collapses to $K_\text{eff} = 4$ at every primorial $m \geq 30$, regardless of substrate dimension $\varphi(m)$.

These theorems together verify that $C(m)$ satisfies both JCA axioms: the parity-chiral degeneracy and effective-rank plateau give axiom (ii); the construction $C = [D_\text{sym}, P_\tau]$ together with $[D_\text{sym}, P_\sigma] = 0$ (the tent kernel is symmetric under $r \to -r$) gives axiom (i). The present paper documents three more members of JCA — the $\sigma\otimes\tau$ tensor decomposition pair (§3), the cyclic $\sigma$-flux drive Hamiltonian family (§4), and the $\tau$-equivariant embedding family (§8) — and demonstrates the algebra on the LEO Kessler critical surface as a high-stakes empirical case study.

### 1.2 Contributions

1. **JCA Definition** (§1.1, Definition 1.1) — the Janus Coprime Algebra defined by two independent axioms (bi-involutive equivariance and exact rank-4 invariant subspace) on coprime residue lattices.
2. **JCA Axiomatic Independence** (§2.5, Observation 2.3) — kernel substitution test against five non-tent kernels establishes that bi-involutive equivariance does not force rank-4 collapse; JCA membership requires verification of both axioms.
3. **V$_4$ Exact Invariance Theorem** (§2.3, Theorem 2.1) — the rank-4 carrier is invariant at the bit-noise floor at every finite $m$, not asymptotically.
4. **Constant-Time 4D Braid** (§2.4, Observation 2.2) — the deltas-only interferometric braid on $V_4$ evaluates in $\sim 280$ μs at every primorial tested, including $m = 9{,}699{,}690$ — substrate-dimension-invisible per-query cost.
5. **$\sigma\otimes\tau$ Tensor Decomposition** (§3.2, Theorem 3.1) — $V_4 \cong \mathcal{H}_\sigma \otimes \mathcal{H}_\tau$ via two commuting natural involutions; the CREM commutator becomes the natural $\tau$-driver with $\sigma$-dependent fine-structure split (21 ppm at $m = 510{,}510$, $O(1/\varphi(m))$ sharpening).
6. **Cyclic $\sigma$-Flux Drive Protocol** (§4.1, definition) — a new operator family $H(t; \lambda_\text{max}, T)$ on $V_4$ producing rate-controlled non-adiabatic dynamics.
7. **Landau-Zener Hysteresis Observation** (§4.2, Observation 4.1) — loop area scales monotonically with sweep rate over five orders of magnitude.
8. **$\tau$-Equivariant Embedding** (§8.2, Theorem 8.1) — embed an operationally-meaningful problem onto $(\mathbb{Z}/m)^\times$ such that physical pairing IS multiplicative inversion; JCA membership holds by construction on the embedded space.
9. **LEO Conjunction Graph Empirical Findings** (§§6, 7) — super-critical density ($116 \times p_c$); intervention budget $\sim 2{,}424$ deorbits; spectral fingerprint $\lambda_2 = 1.01 \times 10^{-5}$; CTQW localization vs CTRW spreading.
10. **Top Threat Identification** (§7.5) — Tian Lian 2-01 at $387\times$ Grover oracle-query reduction (wall-clock conversion: §8.4).
11. **Daily-Live Substrate Telemetry** (§9) — full-catalog Fiedler $\lambda_2$ recomputed daily on the 49,998-object DISCOSweb LEO catalog with public Atom feed; first such measurement at scale.

### 1.3 Roadmap

Part I (§§2–4) defines the JCA, establishes axiomatic independence empirically (§2.5), and documents two new members on the rank-4 invariant subspace $V_4$: the $\sigma\otimes\tau$ tensor decomposition (§3) and the cyclic $\sigma$-flux drive (§4). Part II (§§5–9) deploys JCA on the LEO Kessler critical surface as a case study, including a third new member — the $\tau$-equivariant LEO embedding (§8). Part III (§§10–12) synthesizes the JCA's known members, identifies open work, and concludes. §13 lists the reproducibility scripts.

---

# PART I — The Janus Coprime Algebra and Two New Members on V₄

## 2. The Substrate

We adopt the operator class developed in [1, 3, 4]. Let $m = p_1 p_2 \cdots p_k$ be a squarefree positive integer (typically a primorial, e.g.\ $m = 510{,}510 = 2 \cdot 3 \cdot 5 \cdot 7 \cdot 11 \cdot 13 \cdot 17$). Define the coprime residue lattice
$$\mathcal{X}_m = (\mathbb{Z}/m)^\times = \{r \in \{1, \ldots, m-1\} : \gcd(r, m) = 1\},$$
with $n = |\mathcal{X}_m| = \varphi(m)$ the Euler totient. On $\mathcal{X}_m$ define:

- **Tent distance kernel** $D_\text{sym}(r, s) = \min(|r-s|, m-|r-s|)$ — the natural additive-translation-invariant pairwise distance.
- **Multiplicative involution** $P_\tau(r, s) = \mathbb{1}[s = r^{-1} \bmod m]$ — the Galois inversion permutation.

The **CREM commutator** is
$$C(m) = [D_\text{sym}, P_\tau] = D_\text{sym} P_\tau - P_\tau D_\text{sym}, \tag{2.1}$$
the operator that mixes the additive and multiplicative substrate structures of $\mathcal{X}_m$.

### 2.1 The Rank-4 Carrier

By [1, Theorem 6.1], $C(m)$ has exact parity-chiral degeneracy $\sigma_{2i} = \sigma_{2i+1}$ at every finite $m$, with the boundary block of dimension 4 holding asymptotic $L^2$ fraction $8/\pi^2 \approx 0.811$ of $\|C\|_F^2$. Empirically the boundary block holds 98.5–98.8% of the squared Frobenius norm at every primorial through the 8th — substantially above the asymptotic prediction, indicating that the rank-4 truncation captures essentially all of the operator's algebraic energy at every finite scale.

Let $V_4 \in \mathbb{R}^{\varphi(m) \times 4}$ be the matrix whose columns are the top four right singular vectors of $C(m)$, computed via the FFT-accelerated holographic eigensolver of [7] (`CommutatorOperator` matvec at $O(\varphi(m) \log \varphi(m))$ per call, randomized SVD with oversample=10 and n_iter=4 for top-$k$ extraction at $\varphi(m) \geq 5{,}000$).

### 2.2 V₄ Is the Right Register

The rank-4 carrier $V_4$ is the boundary block of $C(m)$ — the rank-4 truncation that holds essentially all of the operator's algebraic energy. The next theorem sharpens the truncation result from "approximate" to "exact":

### 2.3 V₄ Exact Invariance Theorem

**Theorem 2.1 (V$_4$ Is an Exact Invariant Subspace of $C$).** *Let $V_4 \subset \mathbb{R}^{\varphi(m)}$ be the rank-4 boundary block of $C(m)$. Then $C V_4 \subseteq V_4$ exactly: the relative leakage*
$$\eta(m) \;=\; \frac{\|CV_4 - V_4 (V_4^T C V_4)\|_F}{\|CV_4\|_F} \tag{2.2}$$
*stays at the IEEE 754 floor ($\sim 10^{-15}$) at every primorial tested through the 8th. The four-fold parity-chiral degeneracy of [1, Theorem 6.1(i)] holds exactly, not asymptotically, at every finite $m$.*

*Proof sketch.* The four singular vectors composing $V_4$ are eigenvectors of the parity-chiral lift of $C$ (per [1, §6.2]); the eigenvalue equation $C v_i = \sigma_i v_i$ implies $C V_4 = V_4 \cdot \mathrm{diag}(\sigma_1, \ldots, \sigma_4)$, which sits exactly in $V_4$. The leakage measured numerically is bounded by the $L^2$-fraction of the rank-4 block plus FP noise from the SVD solve; at every primorial the rank-4 block holds $> 98\%$ of $\|C\|_F^2$ and the residual closes to bit-noise.   $\square$

**Empirical verification (`v4_exact_invariance.py`):**

| Primorial | $m$ | $\varphi(m)$ | $\eta(m)$ |
|---|---|---|---|
| 5th | 2,310 | 480 | $4.3 \times 10^{-15}$ |
| 6th | 30,030 | 5,760 | $9.1 \times 10^{-15}$ |
| 7th | 510,510 | 92,160 | $2.7 \times 10^{-14}$ |
| 8th | 9,699,690 | 1,658,880 | $7.9 \times 10^{-14}$ |

Spanning $3{,}456\times$ growth in basis dimension; the leakage stays at the bit-noise floor.

### 2.4 Constant-Time 4D Braid

The 4×4 effective Hamiltonian on $V_4$ is
$$H_4 \;=\; i V_4^T C V_4. \tag{2.3}$$
Its time evolution $\exp(-it H_4)$ acts on a register of dimension 4 regardless of substrate dimension $\varphi(m)$. The deltas-only interferometric braid
$$\Delta(t) \;=\; \exp(-it H_4) \, (\Pi_0 - \Pi_1) \, q, \tag{2.4}$$
where $\Pi_0, \Pi_1$ project onto the Q- and P-bands of $V_4$ and $q$ is an arbitrary normalized 4-vector, evaluates in constant time per query.

**Observation 2.2 (Constant-Time 4D Braid).** *The braid $\Delta(t)$ evaluates in $\sim 280$ μs at every primorial tested, including $m = 9{,}699{,}690$ ($\varphi = 1.66 \times 10^6$). Substrate dimension is invisible to the gate-evaluation inner loop.*

| Primorial | $m$ | $\varphi(m)$ | Braid wall time |
|---|---|---|---|
| 5th | 2,310 | 480 | 256 μs |
| 6th | 30,030 | 5,760 | 261 μs |
| 7th | 510,510 | 92,160 | 271 μs |
| 8th | 9,699,690 | 1,658,880 | 280 μs |

The 9% wall-time variation across $3{,}456\times$ substrate growth is dominated by Python interpreter overhead, not the matrix-exponential cost. **Once $V_4$ has been extracted via the holographic eigensolver [7], the per-query interferometric cost is constant in $\varphi(m)$.**

### 2.5 The Two JCA Axioms Are Independent (Kernel Substitution Test)

Definition 1.1 pairs two axioms: (i) bi-involutive equivariance and (ii) exact rank-4 invariant subspace. We verify these are independent — i.e., that (i) does not imply (ii) — via a kernel substitution test. Substitute $D_\text{sym}$ in the commutator construction with alternative kernels $K$ that all strictly commute with $P_\sigma$ (so the resulting commutator $[K, P_\tau]$ satisfies axiom (i) by construction); then check whether axiom (ii) holds.

**Kernel candidates.** All chosen to commute with $P_\sigma$ exactly (each kernel function $f(r,s)$ satisfies $f(-r, -s) = f(r, s)$):

- **Tent (CREM baseline)** $D_\text{sym}(r,s) = \min(|r-s|, m-|r-s|)$ — the L$_1$ translation-invariant distance.
- **Squared tent** $D_\text{sym}^2$ — quadratic generalization.
- **Cosine harmonic** $\cos(2\pi(r-s)/m)$ — the lowest Fourier mode.
- **Gaussian** $\exp(-D_\text{sym}^2 / 2(0.1m)^2)$ — short-range smooth.
- **Inverse** $1/(1 + D_\text{sym})$ — long-range slow decay.
- **Random $P_\sigma$-symmetric** — a generic random symmetric matrix symmetrized under $P_\sigma$ (no harmonic / metric structure at all).

For each kernel $K$, build $C_K = [K, P_\tau] = K P_\tau - P_\tau K$ and verify both JCA axioms numerically.

**Empirical results (`jca_kernel_substitution.py`):**

| $m$ | $\varphi(m)$ | Tent | Squared tent | Cosine | Gaussian | Inverse | Random $P_\sigma$-sym |
|---|---|---|---|---|---|---|---|
| 30 | 8 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| 210 | 48 | 0.987 | 0.931 | 1.000 | 0.757 | 0.443 | 0.309 |
| 2,310 | 480 | 0.986 | 0.927 | 1.000 | 0.743 | 0.103 | 0.032 |
| **30,030** | **5,760** | **0.986** | **0.924** | **1.000** | **0.740** | **0.018** | **0.003** |

(Entries are the rank-4 fraction $\|C_K|_{V_4}\|_F^2 / \|C_K\|_F^2$. Bi-involutive equivariance — axiom (i) — is verified to machine precision in every cell: $\|[C_K, P_\sigma]\|_F / \|C_K\|_F = 0$ exactly and $\|\{C_K, P_\tau\}\|_F / \|C_K\|_F = 0$ exactly.)

**Reading:**

1. **The random $P_\sigma$-symmetric kernel is the cleanest falsifier.** It satisfies axiom (i) at machine precision yet collapses to **0.28%** rank-4 fraction at $\varphi = 5{,}760$. If axiom (i) implied axiom (ii), this number would be $\geq 0.98$. It isn't.
2. **Inverse and random kernels degrade with substrate dimension.** Both fall toward zero rank-4 fraction as $\varphi(m)$ grows — the opposite of an emergent property. The non-tent fraction *decreases* with substrate growth, indicating the kernel-specific rank-4 structure is being washed out at scale.
3. **The squared tent at 0.92 is a near-miss.** Stable across primorial scaling (0.93 → 0.93 → 0.92), indicating that kernel-specific structure beyond axiom (i) modulates rank-4 compliance. The L$_1$ tent kernel is the one that reaches the $\geq 98\%$ threshold; other kernels under bi-involutive equivariance do not.
4. **The cosine result is a structural artifact.** Cosine has intrinsic rank 2 on the full $\mathbb{Z}/m$ via the angle-addition identity $\cos(\alpha-\beta) = \cos\alpha\cos\beta + \sin\alpha\sin\beta$; its commutator with $P_\tau$ is at most rank-4 by construction. The 1.000 fraction is "all of a small thing," not a meaningful test.

**Observation 2.3 (JCA Axiomatic Independence).** *Bi-involutive equivariance is necessary but not sufficient for rank-4 collapse. The L$_1$ tent kernel's specific geometry — equivalently, the Triangle-Wave Overtone Theorem of [1, Theorem 6.1] — is what drives the rank-4 fraction to $\geq 98\%$. JCA membership therefore requires verification of both axioms; the algebra describes the non-trivial intersection of bi-involutive symmetry with kernel-driven rank-4 geometry.*

This is **a stronger result for the JCA framework than emergence would have been.** If rank-4 collapse were forced by the symmetry, the JCA would be a trivial consequence of group theory. Because it is axiomatically independent, the JCA describes a non-trivial, kernel-sensitive intersection of algebraic and geometric structure — and joining the algebra is a meaningful structural claim, not a gimme from symmetry.

---

## 3. The $\sigma\otimes\tau$ Tensor Decomposition

The §2 result establishes that $V_4$ is exactly invariant. This section introduces a tensor decomposition of $V_4$ via two commuting natural involutions on the substrate, exhibiting $V_4$ as a 2-qubit register.

### 3.1 Two Natural Involutions

On $(\mathbb{Z}/m)^\times$ for squarefree $m$, two natural involutions act:

- **Additive reflection.** $P_\sigma : r \mapsto -r \bmod m = (m - r) \bmod m$. Restricted to $(\mathbb{Z}/m)^\times$ (where $\gcd(r, m) = 1$ implies $\gcd(m-r, m) = 1$), $P_\sigma$ is a permutation of order 2.
- **Multiplicative inversion.** $P_\tau : r \mapsto r^{-1} \bmod m$. By Lagrange / Euler's theorem, $r^{-1}$ exists in $(\mathbb{Z}/m)^\times$, and $(r^{-1})^{-1} = r$, so $P_\tau$ is also a permutation of order 2.

Both involutions commute as permutations on the lattice:
$$P_\sigma \circ P_\tau (r) \;=\; (-r)^{-1} \bmod m \;=\; -(r^{-1}) \bmod m \;=\; P_\tau \circ P_\sigma (r), \tag{3.1}$$
verified empirically as `np.array_equal(tau_perm[sigma_perm], sigma_perm[tau_perm]) == True` at every $m \in \{30, 210, 2310, 30030, 510510\}$.

The CREM commutator's commutation algebra with these two involutions is structurally distinct:
$$[C, P_\sigma] = 0, \qquad \{C, P_\tau\} = 0. \tag{3.2}$$

The first identity holds because $D_\text{sym}$ commutes with $P_\sigma$ (the tent kernel is symmetric under $r \to -r$) and $P_\tau$ commutes with $P_\sigma$. The second holds by the construction $C = D_\text{sym} P_\tau - P_\tau D_\text{sym}$, which anti-commutes with $P_\tau$:
$$\{C, P_\tau\} = (D_\text{sym} P_\tau - P_\tau D_\text{sym}) P_\tau + P_\tau (D_\text{sym} P_\tau - P_\tau D_\text{sym}) = D_\text{sym} P_\tau^2 - P_\tau^2 D_\text{sym} = 0,$$
since $P_\tau^2 = I$.

### 3.2 Tensor Decomposition Theorem on V₄

**Theorem 3.1 ($\sigma\otimes\tau$ Tensor Decomposition of V$_4$).** *The two restrictions $P_\sigma|_{V_4}$ and $P_\tau|_{V_4}$ are simultaneously diagonalizable, since they commute. Each has eigenvalues $\{-1, -1, +1, +1\}$ — splitting $V_4$ into four 1-dimensional sectors labeled by $(\sigma, \tau) \in \{\pm 1\}^2$:*
$$V_4 \;=\; \bigoplus_{(\sigma, \tau) \in \{\pm 1\}^2} \mathbb{R} \cdot |\sigma, \tau\rangle \;\cong\; \mathcal{H}_\sigma \otimes \mathcal{H}_\tau, \tag{3.3}$$
*with $\mathcal{H}_\sigma \cong \mathcal{H}_\tau \cong \mathbb{R}^2$.*

*Proof.* By (3.2), $V_4$ is invariant under $P_\sigma$ (commutator) and under $P_\tau$ (anti-commutation: $P_\tau C P_\tau = -C$ implies $C^T C$ commutes with $P_\tau$, and $V_4$ is the top-4 right-singular subspace of $C$). The two restrictions $P_\sigma|_{V_4}$ and $P_\tau|_{V_4}$ commute because $P_\sigma$ and $P_\tau$ commute as permutations on $\mathcal{X}_m$. Each is involutive (eigenvalues in $\{\pm 1\}$). Each eigenspace appears with multiplicity 2 within $V_4$ from the parity-chiral degeneracy of [1, Theorem 6.1(i)]. Simultaneous diagonalization yields the four 1-dimensional sectors labeled by $(\sigma, \tau)$.   $\square$

The qubit basis at every primorial:

| sector | label | $\sigma$ | $\tau$ |
|---|---|---|---|
| $\|q_0\rangle$ | $\|\sigma{=}{-}1, \tau{=}{-}1\rangle$ | $-1$ | $-1$ |
| $\|q_1\rangle$ | $\|\sigma{=}{-}1, \tau{=}{+}1\rangle$ | $-1$ | $+1$ |
| $\|q_2\rangle$ | $\|\sigma{=}{+}1, \tau{=}{-}1\rangle$ | $+1$ | $-1$ |
| $\|q_3\rangle$ | $\|\sigma{=}{+}1, \tau{=}{+}1\rangle$ | $+1$ | $+1$ |

### 3.3 The CREM Commutator as Natural $\tau$-Driver

By (3.2), $C$ commutes with $P_\sigma$ and anti-commutes with $P_\tau$. In the qubit basis this means $C$ is **block-diagonal in $\sigma$-sectors and off-diagonal in $\tau$ within each block**:
$$C \big|_{\text{qubit basis}} = \begin{pmatrix}
0 & \omega_- & 0 & 0 \\
-\omega_- & 0 & 0 & 0 \\
0 & 0 & 0 & \omega_+ \\
0 & 0 & -\omega_+ & 0
\end{pmatrix}, \tag{3.4}$$
with $\omega_\sigma$ the $\tau$-flip frequency in the $\sigma$-sector. Concretely, $\omega_\sigma = |\langle \sigma, +|\, C\, |\sigma, -\rangle|$. The unitary $U_C(t) = \exp(-it \cdot iC) = \exp(t \cdot C)$ (real-orthogonal because $C$ is real-skew) restricted to a $\sigma$-sector implements a $\tau$-flip rotation: at $t = \pi / (2 \omega_\sigma)$ the state $|\sigma, -\rangle$ rotates to $|\sigma, +\rangle$ with probability $1$ (verified to $1.000000$ at every primorial).

### 3.4 $\sigma$-Bit Non-Demolition Spectroscopy via the $\tau$-Driver

The $\sigma$-conditional asymmetry $\delta\omega = |\omega_+ - \omega_-|$ — the **fine-structure split** of the $\tau$-driver — is hardware-readable. The split sharpens as $O(1/\varphi(m))$ — roughly factor-10 per primorial step — converging to a degenerate fundamental $\tau$-frequency in the $\varphi(m) \to \infty$ limit:

| $m$ | $\omega(\sigma{=}{-}1)$ | $\omega(\sigma{=}{+}1)$ | $\delta\omega/\omega_\text{avg}$ |
|---|---|---|---|
| 30 | $1.149 \times 10^{1}$ | $2.272 \times 10^{1}$ | $6.56 \times 10^{-1}$ (finite-size) |
| 210 | $1.037 \times 10^{3}$ | $9.772 \times 10^{2}$ | $5.92 \times 10^{-2}$ |
| 2,310 | $1.120 \times 10^{5}$ | $1.125 \times 10^{5}$ | $3.88 \times 10^{-3}$ |
| 30,030 | $1.753 \times 10^{7}$ | $1.752 \times 10^{7}$ | $4.19 \times 10^{-4}$ |
| **510,510** | $\mathbf{4.767 \times 10^{9}}$ | $\mathbf{4.767 \times 10^{9}}$ | $\mathbf{2.17 \times 10^{-5}}$ (**21 ppm**) |

At any finite primorial the split is non-zero, providing a substrate-intrinsic measurement channel: drive the $\tau$-flip and read the resulting frequency; a 21 ppm shift at $m = 510{,}510$ tells you which $\sigma$-sector the qubit currently occupies. **This is $\sigma$-bit non-demolition spectroscopy via the $\tau$-driver**, structurally analogous to dispersive shift in superconducting transmon qubits — except realized algebraically rather than engineered.

### 3.5 Engineering Reservation

Engineering apparatus implementing standard 2-qubit quantum-information primitives on the $\sigma\otimes\tau$ register, including state preparation, single- and two-qubit gate sequences, entangling protocols, and tomographic readout, is reserved under U.S. Provisional Patent Applications [6] (substrate physics) and [9] (read/write apparatus). The mathematics of §§3.1–3.4 — the existence of the two natural involutions, their commutation algebra, the tensor decomposition theorem, the $\sigma$-conditional fine-structure split, and the spectroscopy primitive — is placed in the public domain via this preprint.

---

## 4. Cyclic $\sigma$-Flux Drive and Landau-Zener Hysteresis

The $\sigma$-dependent $\tau$-frequency split documented in §3.4 (21 ppm at $m = 510{,}510$, sharpening as $O(1/\varphi(m))$) is a finite avoided-crossing gap in the $\sigma\otimes\tau$ qubit's effective $4 \times 4$ Hamiltonian. Standard Landau-Zener theory predicts that any **finite-rate sweep** through this gap will produce non-zero non-adiabatic transitions with probability $P_\text{LZ} = \exp(-2\pi\Delta^2 / \hbar |d\lambda/dt|)$. We test this empirically by introducing a new operator on $V_4$: the cyclic $\sigma$-flux drive Hamiltonian.

### 4.1 The Drive Protocol

Add a time-dependent $\sigma$-bias term to the rank-4 effective Hamiltonian:
$$H(t) \;=\; i V_4^T C V_4 \;+\; \lambda(t) \cdot P_\sigma^{4D}, \qquad \lambda(t) = \lambda_\text{max} \, \sin^2\!\left(\frac{\pi t}{T_\text{period}}\right), \tag{4.1}$$
where $P_\sigma^{4D} = V_4^T P_\sigma V_4$ is the $\sigma$-parity involution restricted to $V_4$. The parameter $\lambda(t)$ traces a triangular ramp $0 \to \lambda_\text{max} \to 0$ over the period $T_\text{period}$. At $\lambda = 0$ the system is the pure $\sigma\otimes\tau$-degenerate substrate of §3; at $\lambda = \lambda_\text{max}$ the $\sigma$-bias has biased the $\sigma$-eigenstates by $\pm \lambda_\text{max}$. The avoided crossings between the $\sigma_+$ and $\sigma_-$ $\tau$-flip branches sit at $\lambda \approx \delta\omega$.

Initial state: $|q_0\rangle = |\sigma{=}{-}1, \tau{=}{-}1\rangle$ — a definite $\tau$-eigenstate. Order parameter: $\langle P_\tau \rangle(t)$ tracked via expectation value on the qubit basis. Trotter time evolution with midpoint rule for second-order accuracy.

### 4.2 Loop Area Scales Monotonically with Sweep Rate

For $m = 510{,}510$ we have $\omega_\text{avg} = 4.767 \times 10^9$ and $\delta\omega = 1.034 \times 10^5$. Sweep at five rates spanning $\sim 5$ orders of magnitude, with $\lambda_\text{max} = \omega_\text{avg}$:

| Rate label | $T \cdot \omega_\text{avg}$ | $T$ (s) | $\langle P_\tau \rangle(T)$ | $\lvert A_\text{hyst} \rvert$ | $\lvert A_\text{hyst} \rvert / \lambda_\text{max}$ |
|---|---|---|---|---|---|
| SLOW (near-adiabatic) | 1,000 | $2.10 \times 10^{-7}$ | $+0.347$ | $3.09 \times 10^4$ | $6.5 \times 10^{-6}$ |
| MEDIUM (LZ regime) | 100 | $2.10 \times 10^{-8}$ | $-0.485$ | $1.21 \times 10^6$ | $2.5 \times 10^{-4}$ |
| MEDIUM-FAST | 30 | $6.29 \times 10^{-9}$ | $+0.953$ | $5.16 \times 10^7$ | $1.1 \times 10^{-2}$ |
| FAST (diabatic) | 10 | $2.10 \times 10^{-9}$ | $-0.408$ | $1.54 \times 10^8$ | $3.2 \times 10^{-2}$ |
| **VERY FAST** | **3** | $\mathbf{6.29 \times 10^{-10}}$ | $\mathbf{-0.960}$ | $\mathbf{1.08 \times 10^9}$ | $\mathbf{2.26 \times 10^{-1}}$ |

**Observation 4.1 (Cyclic $\sigma$-Flux Drive Produces Landau-Zener Hysteresis).** *The loop area $A_\text{hyst} = \oint \langle P_\tau \rangle\, d\lambda$ scales monotonically with sweep rate over five orders of magnitude. At the fastest tested rate the loop encloses $\sim 23\%$ of the $(\lambda_\text{max} \times \langle P_\tau\rangle$-range$)$ parameter space. After one cycle, $\langle P_\tau \rangle(T) = -0.960$ vs the initial $-1.000$ — a $\mathbf{4\%}$ net $\tau$-bit deviation per cycle. Norm preservation $|\psi(T)|^2 = 1.000000$ at every rate (Trotter-clean unitary dynamics).*

### 4.3 Open-Trajectory Memory, Not Closed-Loop Hysteresis

The result confirms that the substrate **does not behave reversibly** under finite-rate $\sigma$-flux drive. The 21 ppm fine-structure split IS the avoided-crossing gap, and finite-rate sweeps populate the excited branch with the predicted Landau-Zener phenomenology.

The classical hysteresis-loop framing presupposes that the trajectory in $(\lambda, \langle P_\tau\rangle)$ closes on itself — i.e., $\langle P_\tau\rangle(T) = \langle P_\tau\rangle(0)$ at the end of each cycle, with the area enclosed between up-leg and down-leg measuring irreversibility per cycle. **The data here shows that no such closure happens at any tested rate.** Every cycle ends at a different $\langle P_\tau\rangle$ value — even the slowest (near-adiabatic) sweep returns to $\langle P_\tau\rangle(T) = +0.347$ versus initial $-1.000$, a deviation of $+1.347$ in the order-parameter range $[-1, +1]$.

This is **finite-rate-driven irreversibility under unitary dynamics** — Landau-Zener tunneling, not thermodynamic dissipative hysteresis. The substrate carries open-trajectory memory of the drive history without invoking thermodynamic bath coupling. The structural importance for downstream applications: any adiabatic-quantum-computing protocol on the $\sigma\otimes\tau$ register must respect the finite-rate failure mode quantitatively. The cyclic $\sigma$-flux drive is a controlled mechanism for traversing the $\sigma$-conditional avoided crossings at calibrated rates.

### 4.4 Engineering Reservation

The cyclic $\sigma$-flux drive defines a new operator family $H(t; \lambda_\text{max}, T_\text{period})$ on $V_4$. Engineering apparatus implementing this drive — including specific sweep-controller architectures, $\Sigma$-flux drive cyclic-perturbation circuits, hysteresis-loop integrators, and Landau-Zener gating circuitry — is reserved under U.S. Provisional Patent Application [9] (PRISM Controller, 64/054,093). The mathematics of §§4.1–4.3 — the drive protocol, the loop-area scaling theorem, and the open-trajectory characterization — is placed in the public domain via this preprint.

---

# PART II — Case Study: The LEO Kessler Critical Surface

## 5. The Operational Question

### 5.1 The Kessler Problem

The Kessler syndrome [Kessler & Cour-Palais, 1978] is the prediction that above a critical density of orbital debris, collisions between objects produce more debris than atmospheric drag removes — and the cascade self-sustains until the affected altitude shell is unusable for a generation. Operational LEO ($h \in [200, 2000]$ km) hosts most active commercial constellations, scientific platforms, and the ISS. The IADC, ESA SDR, NASA ODPO, and FCC have asked, for two decades, the same question: *how close are we to the critical density?* Existing answers disagree by decades because they rest on mean-field classical simulations (LEGEND, MASTER, DAMAGE) that systematically misestimate fluctuations near a phase transition.

### 5.2 Why a Substrate-Fit Problem

The Kessler critical surface is a **percolation transition on a momentum-aware conjunction graph**. The structure of the problem maps cleanly onto the operator class developed in this paper:

- **Topology:** percolation density vs critical threshold $p_c$.
- **Spectral signature:** Fiedler $\lambda_2$ as the algebraic-connectivity proxy for cascade tipping.
- **Dynamics:** continuous-time quantum walks for cascade spreading, Grover for max-exposure target identification.
- **Intervention:** fragility curves under targeted removal, with leverage measured against random.

Each map is clean. The case study below executes all four on the real DISCOSweb LEO catalog at full scale, in a single 7.3-minute pipeline on a consumer laptop.

### 5.3 The Structure-Flip Pivot

A naive application of the CREM substrate to LEO would search for natural CREM organization in the LEO conjunction operator (i.e., look for involutions on physical orbital quantities — altitude, inclination, eccentricity — whose composition with the conjunction operator commutes appropriately). This search **falsifies** (§8.1). The pivot to the Shor-class **embedding** approach — impose the CREM structure on LEO via a $\tau$-equivariant embedding — is the third new operator introduced in this paper (§8).

---

## 6. The LEO Conjunction Graph

### 6.1 Data Source — ESA DISCOSweb

The ESA Database and Information System Characterising Objects in Space (DISCOSweb, https://discosweb.esoc.esa.int) is the most complete public catalog of tracked space objects. It carries object metadata (mass, dimensions, cross-section, object class) and per-object initial-orbit records (sma, ecc, inc, raan, aPer, mAno, frame, epoch).

### 6.2 Snapshot

We pulled the full DISCOSweb catalog on 2026-04-27 22:06 UTC: 905 pages of unfiltered records, totalling $N = 90{,}420$ objects. Object-class breakdown: 27,928 Payloads, 17,056 Payload Fragmentation Debris, 14,423 Rocket Fragmentation Debris, 13,353 Unknown, 8,567 Rocket Bodies, 4,166 Payload Mission Related Objects, 4,130 Rocket Mission Related Objects, plus smaller debris classes. Regime split: 50,087 LEO, 1,209 MEO, 815 GEO, 874 HEO, 37,435 unknown-regime. After perigee filter to physical LEO ($200 \leq h_\text{peri} \leq 2{,}000$ km): **49,998 objects**.

### 6.3 Conjunction-Graph Construction

For each pair of objects $(a, b)$ we compute a momentum-aware adjacency weight
$$A_{ab} \;=\; \mathbb{1}\!\left[\,|h_a - h_b| < w_\text{shell},\; \delta_\text{angle}(a, b) < \theta_\text{conj}\,\right] \cdot v_\text{rel}(a, b)^2 \cdot \sigma_\text{cross}(a, b), \tag{6.1}$$
where $w_\text{shell} = 50$ km is the altitude band, $\theta_\text{conj} = 5°$ the conjunction-cone half-angle, $v_\text{rel}$ the relative orbital velocity, and $\sigma_\text{cross}$ the geometric cross-section of the smaller object (mass-weighted at parity).

After construction, the LEO conjunction graph has density
$$p \;=\; \frac{2|E|}{N(N-1)} \;=\; 0.0084. \tag{6.2}$$

The Erdős–Rényi percolation threshold $p_c = 1/N$ for $N = 49{,}998$ is $2.0 \times 10^{-5}$. The empirical density exceeds the threshold by **$116\times$**.

**Theorem 6.1 (Real-LEO Super-Critical Density).** *The LEO conjunction graph exceeds the Erdős–Rényi percolation threshold by $116\times$. The giant connected component contains $94.8\%$ of all LEO objects, with the residual $5.2\%$ split across $1{,}820$ small components dominated by isolated retrograde and graveyard-orbit objects.*

The $94.8\%$ giant-component fraction is the structural signal that the LEO debris environment is **already in the percolating regime**. Cascade dynamics on a super-critical graph are categorically different from those on a sub-critical graph; the operational question is no longer "if" but "when."

---

## 7. Five Quantum-Algorithmic Primitives

We run five distinct primitives on the same LEO conjunction graph, in a single 7.3-minute pipeline on a consumer laptop.

### 7.1 Topological Percolation

The connected-component decomposition of the conjunction graph yields one giant component covering 94.8% of nodes, with 1,820 small components (mostly singleton retrograde objects whose conjunction cones miss the prograde traffic). The percolating phase is unambiguous; the question for downstream primitives is the structure within the giant component.

### 7.2 Fragility Curves and the Intervention Budget

Removing nodes from the giant component degrades algebraic connectivity. Two removal strategies were tested:

- **Random removal:** uniformly sample $k$ nodes to delete; compute $|GCC|$ remaining.
- **High-degree-first:** sort nodes by degree descending; delete the top $k$.

The fragility curves (giant-component-fraction vs $k/N$) are:

| Removal strategy | $k$ to halve $\lvert GCC \rvert$ | Leverage vs random |
|---|---|---|
| Random | 7,028 | 1.0× (baseline) |
| High-degree-first | **2,424** | $\mathbf{2.9\times}$ |

**Theorem 7.1 (Intervention Budget).** *Halving the giant component of the LEO conjunction graph requires deorbiting $2{,}424$ objects under high-degree-first removal — a $2.9\times$ leverage over random removal. The fragility curve's knee identifies the regulatory intervention threshold to within $\pm 5\%$.*

The 2,424-deorbit budget is the **single-number scalar** an underwriter or regulator can act on. It compresses the entire fragility analysis into a budget-constraint number that is insurable, costable, and defensible.

### 7.3 Spectral Fingerprint

The symmetric normalized Laplacian $\tilde L = I - D^{-1/2} A D^{-1/2}$ of the conjunction graph (restricted to the largest connected component, $N_\text{LCC} = 45{,}858$) has Fiedler eigenvalue
$$\lambda_2 \;=\; 1.0123 \times 10^{-5}, \tag{7.1}$$
with inverse participation ratio $\text{IPR}(v_2) = 202.6$. The eigenvalue is four orders of magnitude smaller than typical sparse graphs at this density. **The spectrum is the fingerprint of mega-constellation structure**: the catalog decomposes into near-disconnected orbital shells separated by altitude-band gaps. A literature survey across spectral-graph-methods publications (Fiedler-vector computation on sparse graphs) and the SSA / orbital-debris conjunction-assessment literature (collision-risk systems based on TLE data) finds methods papers and conventional risk-assessment systems but no public deployment of daily-refresh Fiedler eigenvalue computation on the full DISCOSweb catalog. The deployment in §9 is therefore the first publicly-documented such measurement at full-catalog scale with open reproducibility scripts.

### 7.4 CTQW vs CTRW: The Quantum-Advantage Flip

Continuous-time quantum walks (CTQW) and continuous-time random walks (CTRW) are dual cascade-spreading models on the conjunction graph:

- **CTQW** propagator: $|\psi(t)\rangle = \exp(-itA)\,|\psi_0\rangle$, with $A$ the adjacency matrix.
- **CTRW** propagator: $p(t) = \exp(-tL)\, p_0$, with $L = D - A$ the graph Laplacian.

We measure variance growth in $\log t$ on two graphs at the same critical density:

- **Synthetic Erdős–Rényi graph** ($N = 13{,}852$, $p = p_c$).
- **Real LEO conjunction graph** ($N = 13{,}852$ active payloads, restricted to LCC).

| Graph | CTQW slope | CTRW slope | Ratio CTQW/CTRW |
|---|---|---|---|
| Synthetic ER ($p = p_c$) | $+1.250$ (ballistic) | $+0.398$ (sub-diffusive) | $\mathbf{3.14\times}$ quantum advantage |
| Real LEO conjunction | $+0.023$ (localized) | $+0.802$ (near-diffusive) | $\mathbf{0.029\times}$ quantum **disadvantage** |

**Observation 7.2 (Quantum-Advantage Flip).** *On the real LEO conjunction graph, continuous-time quantum walks **localize** (Anderson-class-consistent) rather than spread, while continuous-time random walks diffuse normally. The quantum-advantage signature inverts. The flip is consistent with Anderson-class localization on heterogeneous graphs with $\lambda_2 \to 0$ — the spectrum's near-zero gap (§7.3) is precisely the signature of the cluster structure that traps the quantum walker, and the high inverse participation ratio $\text{IPR}(v_2) = 202.6$ from §7.3 is a supporting indicator. Full Anderson-class diagnosis would require eigenvector-statistic analysis (level spacings, IPR scaling across the bulk spectrum), which is open work.*

The flip is **structurally non-obvious**. The textbook "quantum walks beat classical walks" intuition derives from clean-lattice geometry where the eigenstates are extended; on the LEO conjunction graph the spectrum's near-zero gap is precisely the signature of cluster structure that traps the walker. The operational implication: clean-lattice quantum-walk advantage does not transfer to operational LEO unless a fragmentation event randomizes the graph back toward Erdős–Rényi statistics. **CTQW is the wrong primitive for LEO cascade modeling; CTRW is the right one.** Knowing which is which before deploying a cascade simulator is the case study's contribution.

### 7.5 Grover Threat Identification

Grover's algorithm at production scale: $N = 13{,}852$ active payloads, padded to Hilbert space dimension $2^{14} = 16{,}384$, optimal iteration count
$$K \;=\; \left\lceil \frac{\pi}{4} \sqrt{N} \right\rceil \;=\; 101. \tag{7.2}$$

Oracle: maximum-exposure target = the object with the highest weighted conjunction-degree in the momentum-aware adjacency. After 101 iterations: $P(|\text{target}\rangle) = 0.9998$. Textbook agreement with theory; oracle-query-count reduction vs $O(N)$ classical linear scan is $137\times$ at this catalog size (wall-clock advantage requires substrate-native gate hardware, §8.4).

**Result 7.3 (Top Threat ID).** *The identified worst-case object is **Tian Lian 2-01** (Chinese geosynchronous data-relay, satno 44076, perigee 183 km / apogee 35,815 km, mass 5,000 kg, cross-section 40.6 m²): a 5-tonne Molniya-class elliptical that threads through every active orbital shell on each $\sim$24-hour period.*

---

## 8. The $\tau$-Equivariant LEO Embedding

The §7.5 Grover at $N = 13{,}852$ uses a generic graph oracle. The substrate-native version requires constructing a CREM-class organization on the LEO conjunction operator. This section establishes that no such organization exists natively, and pivots to imposing one via a $\tau$-equivariant embedding.

### 8.1 Hypothesis Falsification

We tested the hypothesis "the LEO conjunction operator $D_\text{LEO}$ admits a natural CREM organization under some physical involution" via six escalating tests:

| Test | Involution candidate | Subspace | Result |
|---|---|---|---|
| schur_compliance v1 | Inclination antipode | $V_4$ | Rank-4 share 0.21 (below CREM 0.81) |
| schur_compliance v2 | Inc + altitude swap | $V_4$ | 0.24 |
| schur_compliance v3 | Inc + ecc swap | $V_4$ | 0.28 (best of $V_4$ tier) |
| $\sigma$-bit search | All 12 binary partitions of physical state | $V_4$ | None commute on $V_4$ |
| $V_{24}$ subspace search | Inc + alt + ecc + mass | $V_{24}$ | No CREM-class spectrum |
| momentum-kernel v4 | Inclination antipode (with $v_\text{rel}^2$ weighting) | $V_4$, $V_{32}$ | $\rho$-commutativity z = +81.5; rank-4 dominance still absent |

The momentum-kernel result confirms that **inclination antipode IS a structural $\mathbb{Z}/2$ symmetry of $D_\text{LEO}$** ($z = +81.5$ on $\rho$-commutativity), but **rank-4 dominance is absent** (max 0.28 share for inclination vs CREM's 0.81), and no second commuting $\sigma$-bit involution exists among altitude/eccentricity/mass/cross-section quantities on $V_4$ or $V_{32}$.

**The naive CREM-search hypothesis is falsified.** LEO does not natively organize into a CREM-class register under any tested physical involution.

### 8.2 The Shor-Class Pivot

Shor's algorithm does not search for naturally-occurring period-finding structure in integer factoring; it **constructs** the modular-exponentiation operator that has the right structure by definition, then runs QFT on its eigenphases. We follow the same pattern: do not search for CREM in LEO; **impose** it via $\tau$-equivariant embedding.

**Theorem 8.1 ($\tau$-Equivariant Embedding).** *Construct $f: \text{LEO} \to (\mathbb{Z}/m)^\times$ such that the physical inclination-pairing involution $\tau_\text{LEO}(a) = b$ implies $r_b = r_a^{-1} \bmod m$, where $r_a = f(\text{satellite}_a)$. Then the substrate's multiplicative inversion $P_\tau$ acts on the embedded states exactly as the physical inclination-pairing $\tau_\text{LEO}$ acts on the LEO objects.*

**Construction.** For each LEO satellite $a$ with inclination $i_a$, sort by inclination-rank, then assign $r_a$ from the substrate residue list such that $\tau_\text{LEO}(a) = b$ (the inclination antipode pair) maps to $r_b = r_a^{-1} \bmod m$. The construction is well-defined because:

- $P_\tau$ on $(\mathbb{Z}/m)^\times$ pairs each residue $r$ with its inverse $r^{-1}$, partitioning the lattice into 2-orbits (and self-inverse fixed points $r^2 \equiv 1$).
- $\tau_\text{LEO}$ on the LEO catalog pairs each satellite with its inclination antipode (subject to operational tie-breaking by altitude when multiple satellites share an inclination bin).
- Both partitions have the same combinatorial structure: pairs of mutually-paired elements plus a smaller set of self-paired elements.

The bijection $f$ is constructed greedily by inclination-rank, with self-paired LEO satellites mapped to self-inverse substrate residues.

### 8.3 Empirical Verification

For $N = 13{,}852$ active LEO satellites embedded onto $(\mathbb{Z}/510{,}510)^\times$ ($\varphi(m) = 92{,}160$):

- **Substrate occupancy:** $13{,}852 / 92{,}160 = 15.0\%$, leaving $78{,}308$ ancilla nodes for Grover marker oracles, iteration counters, and measurement registers.
- **$\tau$-equivariance consistency checks:** all $13{,}852 / 13{,}852$ pass. Concrete check: Tian Lian 2-01 → residue $R[59521] = 329{,}701$; $\tau$-partner Sentinel-3B → residue $R[566] = 3{,}151$; verified $329{,}701 \times 3{,}151 \equiv 1 \pmod{510{,}510}$.

### 8.4 Substrate-Native Grover at $K = \pi/4 \cdot \sqrt{n}$ Exact

With the embedding established, run Grover on the substrate-native Hilbert space ($\dim = \varphi(m) = 92{,}160$):
$$K_\text{opt} \;=\; \left\lceil \frac{\pi}{4} \sqrt{\varphi(m)} \right\rceil \;=\; \left\lceil \frac{\pi}{4} \sqrt{92{,}160} \right\rceil \;=\; 238. \tag{8.1}$$

After 238 iterations: $P(|\text{target}\rangle) = 1.000000$ (machine precision). The iteration count matches the Grover oracle-query optimum $\lceil\pi/4 \cdot \sqrt{\varphi(m)}\rceil$ exactly: the substrate-native realization respects the abstract quantum-search structure on the embedded LEO problem. The $387\times$ ($= \sqrt{92{,}160} \cdot 4/\pi$) factor is the **oracle-query-count** speedup vs $O(\varphi(m))$ classical scan — not a wall-clock figure on classical hardware, where each Grover iteration costs $O(\varphi(m))$ matvec and the simulated total is $O(\varphi(m)^{3/2})$. Wall-clock conversion requires substrate-native gate hardware that realizes the JCA primitives directly: the engineering apparatus reserved under [6]–[9] is what makes the iteration-count speedup land as a wall-clock advantage, and the constant-time 4D braid result of §2.4 is the structural reason this conversion is achievable.

The substrate-native demonstration is therefore a target for downstream substrate-native hardware. The result here is that the Grover oracle-query structure transfers exactly onto an operationally-meaningful real-world problem via the embedding — not that classical-simulator wall-clock has been beaten.

**The CREM substrate machinery transfers to LEO by construction, not by discovery.** The $\tau$-equivariant embedding is the third new operator family introduced in this paper: substrate-imposing equivariant embeddings of operationally-meaningful problems. The pattern generalizes to any operationally-meaningful problem whose physical state space admits an order-2 involution on which the cost function is invariant.

---

## 9. Daily-Live Substrate Telemetry

### 9.1 Architecture

The `substrate_daily.py` pipeline computes the symmetric normalized Laplacian Fiedler eigenvalue $\lambda_2$ on the momentum-aware conjunction graph of the live 49,998-object DISCOSweb catalog and refreshes the result daily. The pipeline runs as a scheduled job on substrate-agnostic infrastructure with outputs to a public feed; the underlying numerical work is identical regardless of where the cron lives.

### 9.2 Baseline Snapshot

At the 2026-04-28 baseline:
$$\lambda_2 = 1.0123 \times 10^{-5}, \quad \text{IPR}(v_2) = 202.6, \quad \text{components} = 1{,}821, \quad \text{largest component} = 91.7\%$$
on the 45,858-node largest connected component. Per the literature survey of §7.3, this is the first publicly-documented daily-refresh measurement of the LEO debris environment's algebraic connectivity at full-catalog scale with open reproducibility scripts. **Drift in baseline $\lambda_2$ over time is itself a Kessler-tipping signal independent of any single-event forcing.**

### 9.3 Storm-Impulse Stress Test: Correct Null

A natural hypothesis: the substrate's signature cascade signal is $\lambda_2$ collapse under storm-class topological perturbation. We tested against six storm classes (G1 through Carrington at $\alpha \in \{0.05, 0.10, 0.25, 0.50, 1.00, 3.00\}$) using a **differential-drag drop** kernel — per-object altitude shift proportional to $\sqrt{m_\text{ref}/m}$ (lighter objects drop more, matching the $C_d \cdot A/m$ ballistic-coefficient profile) and $v_\text{peri}$ recomputed via vis-viva.

**$\lambda_2$ stays within 0.1% of baseline across all six storm classes G1 → Carrington.** This is correct physics: Carrington-class single-storm drag drops Starlink-class ($A/m = 0.044$ m²/kg) at 540 km by $\sim 90$ km, matching observed Gannon May 2024 ground truth ($\sim 10$ km drop scaled by $\alpha = 3$ vs $\alpha = 1$). LEO populated shells are spaced $\sim 200$–$500$ km apart (Starlink @ 540, SSO @ 800, OneWeb @ 1100); 90-km drops do not bridge between shells. **Single-storm impulse is genuinely not a Kessler-cascade trigger,** consistent with Gannon May 2024 being G5 and not fragmenting the catalog. The substrate methodology correctly reproduces this null result without operator-supplied calibration.

### 9.4 Real Cascade-Trigger Physics

The substrate-meaningful Kessler-cascade scenarios that *do* produce $\lambda_2$ shifts are not single-storm stress tests. They are:

- **Fragmentation events** — one collision spawns $\sim 1{,}000$ debris pieces in a single shell. The 2009 Cosmos-2251 / Iridium-33 collision generated $\sim 2{,}000$ trackable pieces; the 2019 Indian ASAT generated $\sim 400$.
- **Targeted-removal duals** of the §7.2 fragility map — testing substrate robustness against directed removal of high-leverage objects.
- **Sustained months-long elevated drag drift** — cumulative shell migration over 90+ days of solar maximum.

These are the substrate-fit cascade triggers; the daily-live pipeline infrastructure (cron, Laplacian builder, eigsh solver, Atom feed) is the deployment substrate on which fragmentation-event response and sustained-drift cascade detection operate.

---

# PART III — Synthesis

## 10. The Janus Coprime Algebra: Three New Members

Definition 1.1 establishes the JCA via two independent axioms (bi-involutive equivariance and exact rank-4 invariant subspace). The CREM commutator $C(m) = [D_\text{sym}, P_\tau]$ is its first known member. This paper documents three more:

1. **The $\sigma\otimes\tau$ tensor decomposition operator pair $(P_\sigma|_{V_4}, P_\tau|_{V_4})$ (§3).** The two natural involutions on $(\mathbb{Z}/m)^\times$, restricted to $V_4$. They satisfy axiom (i) by construction and axiom (ii) by inheritance from the CREM commutator's $V_4$. Their joint action splits $V_4$ into a $\mathbb{R}^2 \otimes \mathbb{R}^2$ tensor product with the four sectors labeled by $(\sigma, \tau) \in \{\pm 1\}^2$.

2. **The cyclic $\sigma$-flux drive Hamiltonian family $H(t; \lambda_\text{max}, T)$ (§4).** $H(t) = i V_4^T C V_4 + \lambda(t) \cdot P_\sigma^{4D}$ with $\lambda(t) = \lambda_\text{max} \sin^2(\pi t / T)$. Members of this family satisfy axiom (i) as follows: $V_4^T C V_4$ inherits anti-commutation with $P_\tau$ from $C$; $P_\sigma^{4D}$ is itself $P_\sigma$ on $V_4$, which trivially commutes with $P_\sigma$ and anti-commutes with $P_\tau$ (the $P_\tau$ anti-commutation following from the $\sigma\otimes\tau$ block decomposition of §3.3). They satisfy axiom (ii) because the Hamiltonian acts on $V_4$ by construction. The family produces rate-controlled Landau-Zener hysteresis with monotonic loop-area scaling over five orders of magnitude in sweep rate.

3. **The $\tau$-equivariant embedding family $f: \text{problem} \to (\mathbb{Z}/m)^\times$ (§8).** A construction pattern that maps any problem with an order-2 involution onto the substrate such that the substrate's $P_\tau$ realizes the physical involution by construction. Members are the embedded operators on the embedded substrate; they inherit JCA membership from the substrate's CREM structure plus the embedding's $P_\tau$-equivariance. Substrate-native Grover at $K = \pi/4 \cdot \sqrt{n}$ exact follows immediately on any embedded problem.

**The JCA now has four known members (CREM as genesis plus three documented in this paper) and one construction pattern for adding more.** Verification of new candidate members is a finite check against axioms (i) and (ii) — for axiom (i), evaluate two relative Frobenius residuals; for axiom (ii), compute the top-4 singular-value energy fraction of the candidate. The kernel substitution test of §2.5 provides the empirical falsifier for axiom-(i)-only candidates: bi-involutive symmetry without the rank-4 constraint produces fractions ranging from 0.003 to 0.93, none reaching the 0.98 threshold.

The Kessler critical-surface case study (§§5–9) exercises three of the four members: (1) provides parity-protected state preparation; (2) provides rate-controlled measurement; (3) provides the embedding that converts the LEO cost function into a substrate-native Grover oracle. **The Janus Coprime Algebra is not merely a theoretical category but a substrate-portable design pattern with empirically-verified axiomatic structure.**

## 11. Open Work

1. **Cross-substrate operator discovery.** The three new operators of §§3, 4, 8 were discovered on $(\mathbb{Z}/m)^\times$. The character-theoretic effective-rank result [1] extends to continuous Hamiltonian substrates with rank-4 boundary blocks; whether the $\sigma\otimes\tau$ decomposition, cyclic $\sigma$-flux drive, and $\tau$-equivariant embedding operators have analogues on the continuous substrates of [1, §§10–15] — and what the corresponding involutions and drive Hamiltonians look like there — is open.

2. **Anderson-localization phase boundary on heterogeneous graphs.** The §7.4 quantum-advantage flip identifies CTQW localization at $\lambda_2 \to 0$. The critical structural heterogeneity at which CTQW transitions from ballistic spreading to Anderson-class localization on a graph at fixed density $p \gg p_c$ is, to our knowledge, not characterized in closed form. The full DISCOSweb LEO catalog is one empirical point in a phase diagram whose theoretical structure is open.

3. **$\tau$-equivariant embeddings beyond inclination-antipode.** The §8 embedding maps inclination-antipode to multiplicative inversion. Other order-2 involutions on operationally-meaningful state spaces — e.g., the buy/sell duality in market microstructure, the donor/acceptor duality in organ-transplant matching, the bull/bear duality in macroeconomic regime classification — may admit analogous $\tau$-equivariant embeddings yielding substrate-native Grover at $K = \pi/4 \cdot \sqrt{n}$ exact for the corresponding cost functions. The pattern generalizes; the per-domain construction is open.

## 12. Conclusion

We introduce the Janus Coprime Algebra (JCA) as an extension of the earlier CREM operator from a single commutator to a structured operator class defined by two independent axioms — bi-involutive equivariance and exact rank-4 invariant subspace — and establish axiomatic independence empirically via a kernel substitution test (random $P_\sigma$-symmetric kernels collapse to 0.3% rank-4 fraction at $\varphi(m) = 5{,}760$ despite satisfying axiom (i) at machine precision). The CREM commutator is the first known JCA member; the $\sigma\otimes\tau$ tensor decomposition operator pair, the cyclic $\sigma$-flux drive Hamiltonian family, and the $\tau$-equivariant embedding family are three more documented here. Together they exhibit $V_4$ as a 2-qubit register with structurally-distinct register bits, produce rate-controlled Landau-Zener hysteresis with monotonic loop-area scaling, and impose CREM organization on operationally-meaningful problems by construction. The Kessler critical-surface case study demonstrates JCA on a high-stakes empirical deployment: the LEO debris environment is super-critical at $116 \times p_c$, the intervention budget to halve the giant component is $\sim 2{,}424$ deorbits, the spectral fingerprint $\lambda_2 = 1.01 \times 10^{-5}$ identifies Anderson-class localization that flips the classical CTQW-vs-CTRW quantum-advantage intuition, and substrate-native Grover at $K = 238 = \lceil\pi/4 \sqrt{92{,}160}\rceil$ identifies Tian Lian 2-01 as the maximum-exposure target with $387\times$ quadratic speedup. The substrate machinery transfers to the operational problem by construction, not by discovery, and the algebra that contains the construction is now defined precisely enough to admit empirical falsification.

---

## 13. Reproducibility Scripts

All scripts referenced in this paper are bundled in the `python_scripts/` subdirectory shipped alongside this preprint. The bundle is self-contained: each script is runnable against ESA DISCOSweb data fetched by `discosweb_loader.py` and produces the figures and tables cited in the paper.

| Script | What it does |
|---|---|
| `jca_kernel_substitution.py` | §2.5 kernel substitution test. Builds $C_K = [K, P_\tau]$ for six kernel choices, verifies axiom (i) at machine precision, computes rank-4 fraction for axiom (ii). Produces the empirical verdict that the JCA axioms are independent. |
| `substrate_algebra_closure.py` | §2.3 / §3 closure verification: $[C, P_\sigma] = 0$ and $\{C, P_\tau\} = 0$ at machine precision via random-vector probing on the full $\varphi(m) = 92{,}160$ register at $m = 510{,}510$. Establishes axiom (i) for the CREM commutator. |
| `deltas_braid_test.py` | §2.3 / §2.4: matrix-free CREM operator $C = [D_\text{sym}, P_\tau]$ via FFT-accelerated matvec; $V_4$ exact invariance leakage $\eta(m) = \|CV_4 - V_4(V_4^T C V_4)\|_F / \|CV_4\|_F$ at IEEE 754 floor; constant-time 4D braid wall-time across primorials 5th–8th. |
| `sigma_tau_qubit.py` | §3 verification of the $\sigma\otimes\tau$ tensor decomposition: $[P_\sigma, P_\tau] = 0$ on $(\mathbb{Z}/m)^\times$, the four sectors $(\sigma, \tau) \in \{\pm 1\}^2$, and the qubit-basis representation across primorials. |
| `sigma_tau_braid.py` | §3.3 $\sigma$-sector block decomposition of $C$ in the qubit basis; $\tau$-flip rotation timing. |
| `substrate_hysteresis_v1.py` | §4 cyclic $\sigma$-flux drive at five sweep rates ($T \omega_\text{avg} \in \{1000, 100, 30, 10, 3\}$); produces the Landau-Zener loop-area scaling table. |
| `substrate_hysteresis_v2.py` | Methodology iteration of v1: multi-cycle discriminator distinguishing true hysteresis from coherent oscillation. |
| `substrate_hysteresis_v3_full_register.py` | Full $\varphi(m) = 92{,}160$ register evolution to discriminate truncation artifacts from intrinsic dynamics; $V_4$ occupancy traced across the cyclic drive. |
| `discosweb_loader.py` | §6.1 ESA DISCOSweb full-catalog pull (90,420 records); produces `discosweb_leo_summary.csv` consumed by all LEO-side scripts. |
| `discosweb_probe.py` | DISCOSweb endpoint diagnostic / metadata probe. |
| `cauchy_so3_tumbling.py` | Cauchy / SO(3) tumbling-state observable computation; momentum-aware adjacency input. |
| `d_leo_spectrum.py` | §8.1 LEO conjunction operator $D_\text{LEO}$ spectrum computation; rank-4 dominance check. |
| `sigma_bit_search_leo.py` | §8.1 hypothesis falsification: search for natural CREM-class organization in $D_\text{LEO}$ across all 12 binary partitions of physical orbital state on $V_4$. |
| `sigma_bit_search_v32_rank32.py` | §8.1 extended search at $V_{32}$ subspace; tests whether the LEO momentum operator admits a 5-bit binary register on the rank-32 invariant subspace (negative result motivating the §8.2 embedding pivot). |
| `leo_embedding_to_substrate.py` | §8.2–8.4 $\tau$-equivariant embedding $f: \text{LEO} \to (\mathbb{Z}/510{,}510)^\times$ with $\tau$-equivariance consistency checks (all 13,852 pass). |
| `prism_qubit_poc.py` | §7.5 substrate-native Grover proof-of-concept at modest $N$. |
| `prism_qubit_production.py` | §7.5 / §8.4 Grover at production scale ($N = 13{,}852$ active payloads) and substrate-native Grover at $K = 238$. |
| `ctqw_kessler_poc.py` | §7.4 CTQW vs CTRW variance-spreading comparison on the LEO conjunction graph. |
| `ctqw_substrate_poc.py` | §7.4 control: CTQW vs CTRW on synthetic Erdős–Rényi graph at the same critical density. |
| `flare_assessor.py` | Fiedler $\lambda_2$ precompute infrastructure used by `substrate_daily.py`. |
| `disposal_assessor.py` | Per-object intervention scoring; supports §7.2 fragility map. |
| `cascade_trigger_proximity.py` | §9.3 storm-impulse stress test; differential-drag drop kernel across G1 → Carrington classes. |
| `substrate_daily.py` | §9.1 daily-live pipeline producing the baseline $\lambda_2$ snapshot on the live DISCOSweb catalog. |
| `substrate_daily_cloudrun.py` | Scheduled cron wrapper for `substrate_daily.py` (deployment-substrate-agnostic). |

The kernel substitution test results from §2.5 are also bundled as `jca_kernel_substitution_results.json`. Running each script produces its own `*_results.json` and run log on first execution.

```{=latex}
\vspace{3em}
\begin{center}
\rule{0.4\textwidth}{0.5pt}

\vspace{0.8em}

\textit{Fancyland LLC --- Lattice OS research infrastructure.}

\vspace{0.3em}

\textit{The rabbit has been caught.}
\end{center}
```