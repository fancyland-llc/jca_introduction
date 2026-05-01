# Introducing the Janus Coprime Algebra

**A Case Study on the LEO Kessler Critical Surface**

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.19955555-blue)](https://doi.org/10.5281/zenodo.19955555)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the preprint and the reproducibility bundle for *Introducing the Janus Coprime Algebra: A Case Study on the LEO Kessler Critical Surface* (Matos, 2026). The paper extends the earlier CREM operator from a single commutator to a structured operator class — the **Janus Coprime Algebra (JCA)** — defined on coprime residue lattices `(Z/m)^×` for squarefree primorial `m` by two independent axioms (bi-involutive equivariance and exact rank-4 invariant subspace), documents three new members of the algebra, and demonstrates the construction on a high-stakes empirical deployment: locating the Kessler critical surface in the modern LEO debris environment.

The JCA's known members:

1. **CREM commutator** *C(m) = [D_sym, P_τ]* — the genesis member, characterized in the foundational stack.
2. **σ⊗τ tensor decomposition** (§3) — *V_4 ≅ H_σ ⊗ H_τ* via two commuting natural involutions; the CREM commutator becomes the natural τ-driver with σ-dependent fine-structure split (21 ppm at *m = 510,510*, *O(1/φ(m))* sharpening).
3. **Cyclic σ-flux drive Hamiltonian family** (§4) — rate-controlled non-adiabatic dynamics on *V_4* with monotonic loop-area scaling over five orders of magnitude in sweep rate.
4. **τ-equivariant embedding family** (§8) — substrate-imposing equivariant embeddings of operationally-meaningful problems onto `(Z/m)^×`; substrate-native Grover at the oracle-query optimum follows by construction.

Five named theorems anchor the work:

1. **JCA Definition** (§1.1) — two independent axioms (bi-involutive equivariance + exact rank-4 invariant subspace) on coprime residue lattices.
2. **JCA Axiomatic Independence** (§2.5) — kernel substitution test against five non-tent kernels establishes that bi-involutive equivariance does not force rank-4 collapse; the L₁ tent kernel's specific geometry is what drives the rank-4 fraction to ≥ 98%. JCA membership requires verification of both axioms.
3. **V₄ Exact Invariance Theorem** (§2.3) — *η(m) = ‖CV_4 − V_4(V_4ᵀCV_4)‖_F / ‖CV_4‖_F* stays at the IEEE 754 floor (~10⁻¹⁵) at every primorial tested through the 8th, spanning 3,456× growth in basis dimension. The four-fold parity-chiral degeneracy holds exactly, not asymptotically.
4. **σ⊗τ Tensor Decomposition Theorem** (§3.2) — the two natural involutions *P_σ* (additive reflection) and *P_τ* (multiplicative inversion) commute as permutations and split *V_4* into four 1-dimensional sectors labeled by *(σ, τ) ∈ {±1}²*.
5. **τ-Equivariant Embedding Theorem** (§8.2) — embed any problem with an order-2 involution onto `(Z/m)^×` such that the substrate's *P_τ* realizes the physical involution by construction; the JCA structure transfers by construction, not by discovery.

Empirical verification spans the rank-4 carrier across primorials 5–8 (substrates *m ∈ {2310, 30030, 510510, 9699690}* with *φ(m) ∈ {480, 5760, 92160, 1658880}*), the full ESA DISCOSweb LEO catalog (49,998 objects after physical-LEO perigee filter), and the LEO embedding (13,852 active payloads embedded into `(Z/510510)^×`).

## Repository Structure

```
paper/
  JCA_INTRODUCTION.pdf     ← Preprint PDF (the artifact archived on Zenodo)
  JCA_INTRODUCTION.md      ← Markdown source
  JCA_INTRODUCTION.tex     ← LaTeX source (pandoc-generated)
scripts/
  jca_kernel_substitution.py        ← §2.5 axiomatic independence test
  substrate_algebra_closure.py      ← §2.3 / §3 closure verification
  deltas_braid_test.py              ← §2.3 / §2.4 V_4 invariance + 4D braid
  sigma_tau_qubit.py                ← §3 σ⊗τ tensor decomposition
  sigma_tau_braid.py                ← §3.3 σ-sector block decomposition
  substrate_hysteresis_v1.py        ← §4 Landau-Zener loop-area scan (5 rates)
  substrate_hysteresis_v2.py        ← §4 multi-cycle discriminator
  substrate_hysteresis_v3_full_register.py ← §4 full-register V_4 occupancy
  discosweb_loader.py               ← §6.1 ESA DISCOSweb catalog pull
  discosweb_probe.py                ← DISCOSweb endpoint diagnostic
  cauchy_so3_tumbling.py            ← Cauchy / SO(3) tumbling-state observable
  d_leo_spectrum.py                 ← §8.1 LEO conjunction operator spectrum
  sigma_bit_search_leo.py           ← §8.1 hypothesis falsification (V_4)
  sigma_bit_search_v32_rank32.py    ← §8.1 extended search (V_32)
  leo_embedding_to_substrate.py     ← §8.2-8.4 τ-equivariant LEO embedding
  prism_qubit_poc.py                ← §7.5 substrate-native Grover proof-of-concept
  prism_qubit_production.py         ← §7.5 / §8.4 Grover at production scale
  ctqw_kessler_poc.py               ← §7.4 CTQW vs CTRW on LEO
  ctqw_substrate_poc.py             ← §7.4 control on synthetic Erdős-Rényi
  flare_assessor.py                 ← Fiedler λ_2 precompute infrastructure
  disposal_assessor.py              ← §7.2 fragility per-object scoring
  cascade_trigger_proximity.py      ← §9.3 storm-impulse stress test
  substrate_daily.py                ← §9.1 daily-live pipeline
  substrate_daily_cloudrun.py       ← Scheduled cron wrapper (deployment-substrate-agnostic)
  jca_kernel_substitution_results.json ← Cached §2.5 results
LICENSE                             ← MIT
```

## Reproducing the Theorems

Requires **Python ≥ 3.10**, **NumPy ≥ 2.0**, **SciPy ≥ 1.11**. The DISCOSweb-side scripts also require `requests` and a free DISCOSweb API token (set as `DISCOSWEB_TOKEN` environment variable; tokens are issued at <https://discosweb.esoc.esa.int/tokens>). No GPU. No proprietary services.

### Part I — Janus Coprime Algebra (substrate-only)

```bash
cd scripts

# §2.5 — axiomatic independence (six kernels, four primorials)
python jca_kernel_substitution.py

# §2.3 / §3 — bi-involutive closure verification
python substrate_algebra_closure.py

# §2.3 / §2.4 — V_4 exact invariance + constant-time 4D braid
python deltas_braid_test.py

# §3 — σ⊗τ tensor decomposition + qubit basis
python sigma_tau_qubit.py
python sigma_tau_braid.py

# §4 — cyclic σ-flux drive, Landau-Zener loop-area scaling
python substrate_hysteresis_v1.py
python substrate_hysteresis_v3_full_register.py
```

Each script writes a `*_results.json` artifact and a console run log.

### Part II — LEO Kessler Critical Surface

```bash
cd scripts

# §6.1 — pull DISCOSweb catalog (requires DISCOSWEB_TOKEN env var)
python discosweb_loader.py --query leo

# §8.1 — natural-CREM hypothesis falsification on D_LEO
python d_leo_spectrum.py
python sigma_bit_search_leo.py
python sigma_bit_search_v32_rank32.py

# §8.2-8.4 — τ-equivariant LEO embedding + substrate-native Grover
python leo_embedding_to_substrate.py
python prism_qubit_production.py

# §7.4 — CTQW vs CTRW quantum-advantage flip
python ctqw_substrate_poc.py    # synthetic Erdős-Rényi control
python ctqw_kessler_poc.py      # real LEO conjunction graph

# §7.2 — fragility curves and intervention budget
python disposal_assessor.py

# §9 — daily-live substrate telemetry
python substrate_daily.py

# §9.3 — storm-impulse stress test (six storm classes)
python cascade_trigger_proximity.py
```

End-to-end Part II wall time: approximately 7 minutes on a Ryzen 9 7950X single core, dominated by the Fiedler eigenvalue solve on the 45,858-node largest connected component.

## Companion Patents

The mathematics in this preprint underlies a **trilogy+1** of US provisional patent applications protecting the corresponding engineering apparatus. The mathematics is published openly under MIT license; the engineering apparatus, system embodiments, gate-sequence implementations, and method claims drawn from this mathematics are reserved.

| Layer | US Provisional | Title | Filed |
|---|---|---|---|
| State storage | **64/031,440** | Fault-Injection-Immune Computational Unit Using Primorial Coprime Residue Topology (the Arithmetic Qubit) | 2026-04-07 |
| Spectral compute | **64/033,689** | Holographic Eigen-Solver Using QM Boundary Projection on Coprime Residue Lattices | 2026-04-08 |
| Stereoscopic measurement | **64/048,617** | Tensegrity Interferometer: Stereoscopic Query Resolution and Manifold-Curvature Measurement on Continuous Hamiltonian Substrates | 2026-04-24 |
| Read/write control | **64/054,093** | Method and Apparatus for Birefringent Control of a Quantum Arithmetic Substrate via Bloch-Flux-Driven P-Band Modulation (the PRISM Controller) | 2026-04-30 |

Implementations practicing the disclosed mathematics in the manner claimed in the patents above require a license from the assignee (Fancyland LLC).

## Citation

```bibtex
@misc{matos2026jca,
  author       = {Matos, Antonio P.},
  title        = {Introducing the Janus Coprime Algebra: A Case Study on the LEO Kessler Critical Surface},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.19955555},
  url          = {https://doi.org/10.5281/zenodo.19955555}
}
```

## Companion Preprints

The foundational mathematical stack:

- **Character-Theoretic Effective Rank** — DOI [10.5281/zenodo.19744573](https://doi.org/10.5281/zenodo.19744573) ([github.com/fancyland-llc/character-rank](https://github.com/fancyland-llc/character-rank))
- **The Arithmetic Black Hole: Softmax Thermodynamics and the Four Eigenvalue Laws of the Prime Gas** — DOI [10.5281/zenodo.19442006](https://doi.org/10.5281/zenodo.19442006)
- **The Unity Clock: Effective Dimensional Collapse of the Addition-Multiplication Commutator on Coprime Lattices** — DOI [10.5281/zenodo.19478727](https://doi.org/10.5281/zenodo.19478727)
- **Active Transport on the Prime Gas: Flat-Band Condensation, the Rabi Phase Transition, and the Arithmetic Qubit** — DOI [10.5281/zenodo.19243258](https://doi.org/10.5281/zenodo.19243258)
- **Universal Two-Prime Formula for the Coprime-Lattice Coupling Constant** — DOI [10.5281/zenodo.19210625](https://doi.org/10.5281/zenodo.19210625)

---

*Fancyland LLC — Lattice OS research infrastructure.*
