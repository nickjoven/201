# Gravity as Synchronization in a Frictional Medium

**N. Joven** &middot; [ORCID: 0009-0008-0679-0812](https://orcid.org/0009-0008-0679-0812) &middot; March 2026

A unified framework proposing that gravity is not a force of attraction but a
phenomenon of **synchronization** through a **frictional medium**. Two tabletop
experiments — a bowed brass tube and metronomes on a swing — motivate the full
theoretical structure, from galaxy rotation curves to black holes to the CMB.

> **[Read the walkthrough &rarr;](https://nickjoven.github.io/201/)**

---

## Core Thesis

| Concept | Standard Picture | This Framework |
|---|---|---|
| Gravity | Curvature of spacetime by mass | Synchronization of oscillators through a frictional medium |
| Metric tensor | Describes curvature | Local friction coefficient of the vacuum |
| Dark matter | Invisible mass | Shadow price (Lagrange multiplier) for failed synchronization |
| MOND transition | Empirical fitting function | Kuramoto coupling threshold of the medium |
| Black hole horizon | Escape velocity = *c* | Duty cycle = 1.0 (permanent stick, no slip) |

---

## Repository Structure

```
201/
├── joven_unifying_framework.md        ← Main paper
├── renzos_rule_derivation.md          ← Formal derivation of Renzo's Rule & inverse
├── kuramoto_einstein_mapping.md       ← Explicit derivative mapping (ADM ↔ Kuramoto)
├── sparc_x/                           ← Python implementation of the framework
│   ├── calculator.py                      Unified API (MOND + Kuramoto paths)
│   ├── kuramoto.py                        Continuum Kuramoto solver on radial manifold
│   ├── lyapunov.py                        Lyapunov functional & stability analysis
│   ├── stribeck.py                        Stribeck friction ↔ MOND mapping
│   ├── mond.py                            MOND interpolating functions (4 families)
│   ├── profiles.py                        SPARC galaxy profile loading & synthetics
│   ├── constants.py                       Physical constants
│   └── tests/                             Test suite
├── data/                              ← 175 SPARC galaxy rotation curves
├── intersections/                     ← Submodule: companion papers & notebooks
│   ├── joven_stick_slip_dark_matter.md    Companion paper (local mechanism, dark matter dual)
│   ├── cone_topology.ipynb                Flat rotation curves from conical geometry
│   ├── stick_slip_galaxy.ipynb            Stick-slip dynamics applied to galaxies
│   ├── stick_slip_lagrangian.ipynb        Lagrangian relaxation derivation
│   ├── 08_square_wave_bifurcation.ipynb   Square-wave bifurcation analysis
│   ├── 09_cylinder_stick_slip.ipynb       Cylindrical stick-slip mode coupling
│   ├── dispersion_rolling.ipynb           Dispersion and rolling analysis
│   ├── cvt/                               Consistency Vector Theory
│   │   ├── synthesis.md                       CVT synthesis document
│   │   ├── laws/noninjectivity.md             Law of Genealogical Non-Injectivity
│   │   └── notebooks/
│   │       ├── 01_sparc_mode_decomposition.ipynb
│   │       ├── 02_gravitational_fundamental.ipynb
│   │       ├── 03_inharmonicity_diagnostic.ipynb
│   │       ├── 04_mode_coupling_cascade.ipynb
│   │       ├── 05_cmb_overtone_comparison.ipynb   CMB peaks as overtone series
│   │       ├── 06_feigenbaum_cascade.ipynb         Period-doubling to black holes
│   │       └── 07_equipartition_uv_cutoff.ipynb
│   └── index.html                         Intersections landing page
├── docs/
│   └── index.html                     ← GitHub Pages walkthrough
├── .ket/                              ← Ket substrate (reserved)
├── LICENSE                            ← CC0 1.0 — No rights reserved
└── README.md                          ← This file
```

---

## Key Results

### Mathematically Established

1. **Stick-slip bifurcation is bidirectional.** Both slow drive and overwhelming
   force produce subharmonics through the same mechanism. Experimentally
   validated by Kawano et al. (2025). See the
   [companion paper](intersections/joven_stick_slip_dark_matter.md) §1–2.

2. **Stribeck friction curve is structurally parallel to MOND interpolation.**
   The velocity-weakening branch maps onto the MOND interpolating function
   &mu;(*x*). See the
   [companion paper](intersections/joven_stick_slip_dark_matter.md) §2.1.

3. **Dark matter halos as Lagrange multipliers.** Constrained optimization
   yields dual variables that reproduce NFW-like profiles when baryonic matter
   is insufficient for the gravitational constraint. See the
   [companion paper](intersections/joven_stick_slip_dark_matter.md) §3–4.

4. **Renzo's Rule from complementary slackness.** Every feature in the baryonic
   luminosity profile is mirrored in the rotation curve as a structural
   consequence of KKT conditions, not an empirical correlation. See the
   [companion paper](intersections/joven_stick_slip_dark_matter.md) §4.3 and
   [formal derivation](renzos_rule_derivation.md).

5. **Inverse Renzo's Rule: no phantom features.** Every localized rotation curve
   feature has a baryonic origin. The Green's function of the Laplacian is a
   smoothing kernel — it cannot sharpen or generate localized structure. Proved
   from Poisson equation alone (no mapping dependence). See the
   [formal derivation](renzos_rule_derivation.md) §3.3. The Lyapunov functional
   implemented in [`sparc_x/lyapunov.py`](sparc_x/lyapunov.py) proves
   *dV*/*dt* &leq; 0 along Kuramoto trajectories, establishing convergence to a
   unique attractor at stable equilibria — upgrading the inverse from "no phantom
   features" to "no phantom features or backgrounds" at physically realized
   fixed points.

6. **Flat rotation curves from conical geometry.** On a conical metric, Bessel
   modes yield *k* &prop; 1/*r* without any dark matter term (&lambda; = 0). See the
   [cone topology notebook](intersections/cone_topology.ipynb).

7. **Kuramoto synchronization threshold is formally identical to the MOND
   transition.** Coupling strength *K<sub>c</sub>* maps to acceleration scale
   *a*<sub>0</sub>. Below threshold: incoherence (Newtonian). Above threshold:
   phase-locking (MOND). See the
   [main paper](joven_unifying_framework.md) §3.3.

8. **Feigenbaum cascade toward black hole singularity.** Period-doubling at the
   universal rate &delta; &approx; 4.669 accumulates at the event horizon. See the
   [cascade notebook](intersections/cvt/notebooks/06_feigenbaum_cascade.ipynb).

### Reasonable Next Hypotheses

Items marked ✓ have been resolved or substantially advanced by the
[`sparc_x`](sparc_x/) implementation. See the detailed status in
[renzos_rule_derivation.md §7](renzos_rule_derivation.md),
[kuramoto_einstein_mapping.md §7](kuramoto_einstein_mapping.md), and
[joven_unifying_framework.md §9](joven_unifying_framework.md).

9. **The metric tensor *is* the local friction coefficient.** The Einstein field
   equations reread as: the medium's impedance (*G*<sub>&mu;&nu;</sub>) is
   determined by its energy content (*T*<sub>&mu;&nu;</sub>). Pressure and
   velocity determine nodal density. This reframes GR with identical
   predictions but different ontology. *Partial progress:* the Stribeck-to-MOND
   correspondence is now [explicit](sparc_x/stribeck.py) and the
   [Kuramoto-Einstein mapping](kuramoto_einstein_mapping.md) establishes
   term-by-term structural correspondence. *Remains:* verifying all numerical
   prefactors (dynamical equivalence).

10. **Space is a phase shift in the medium's duty cycle.** The ADM lapse
   function *N* and shift vector *N<sup>i</sup>* are duty cycle and phase
   offset of local oscillation. *Status:* the ADM formalism already contains
   this structure; the claim is that the physical interpretation is literal,
   not merely formal.

11. **The Universal Rosin — a measurable vacuum impedance function.** The
    vacuum's Stribeck curve sets *a*<sub>0</sub> at galactic scales, *I<sub>c</sub>*
    at Josephson junctions, and rate-and-state friction at fault zones. *Needs:*
    independent measurement of the vacuum Stribeck curve, possibly via the
    relationship between &Lambda; and *a*<sub>0</sub>.

12. **CMB acoustic peaks as primordial synchronization.** Peaks mark scales
    where oscillator synchronization completed at decoupling; the Sachs-Wolfe
    plateau marks scales too large to synchronize. See the
    [CMB notebook](intersections/cvt/notebooks/05_cmb_overtone_comparison.ipynb).
    *Needs:* sub-percent accuracy match to Planck data.

13. **The Bullet Cluster as synchronization decoupling.** During violent
    mergers, synchronization tracks the collisionless component while gas
    decouples. Plausible within the framework but undemonstrated.

---

## Open Questions in Other Frameworks This Resolves

| Open Question | Framework | Resolution |
|---|---|---|
| Why does MOND work so well empirically? | &Lambda;CDM | MOND transition *is* the Kuramoto synchronization threshold — it has structural origin, not coincidental fit |
| Why does dark matter track baryonic features (Renzo's Rule)? | &Lambda;CDM | Complementary slackness: dual variable is structurally coupled to primal variable at every radius. Inverse also holds: at the Kuramoto-Einstein fixed point, &rho;<sub>dark</sub> is a functional of &rho;<sub>bary</sub> |
| What physical mechanism produces *a*<sub>0</sub>? | MOND | *a*<sub>0</sub> is the impedance-matching point of the vacuum's Stribeck curve — the scale where medium transitions from stiff to compliant |
| Why do galaxy rotation curves flatten? | Both | On conical geometry, Bessel modes yield *k* &prop; 1/*r* automatically; or equivalently, synchronization coupling produces flat coherence beyond the threshold radius |
| Why are QPO frequency ratios approximately mass-independent? | GR | Synchronization ratios are set by medium nodal geometry, not by gravitational source mass |
| What is dark matter? | &Lambda;CDM | Shadow price (Lagrange multiplier) the gravitational field pays when baryonic matter is insufficient for synchronization — it is constraint structure, not substance |
| What happens at a singularity? | GR | Feigenbaum cascade accumulation point: the medium transitions from oscillation to permanent stick — a phase transition, not a breakdown of physics |

---

## Companion Repository

The [intersections](https://github.com/nickjoven/intersections) repository
contains the companion paper, computational notebooks, and the Consistency
Vector Theory framework referenced throughout. It is included here as a git
submodule. To populate it after cloning:

```bash
git submodule update --init --recursive
```

For additional context, see also
[proslambenomenos](https://github.com/nickjoven/proslambenomenos) — an
independent reference framework with complementary perspective. This
repository (201) is fully self-contained and does not depend on it.

---

## License

[CC0 1.0 Universal](LICENSE) — No rights reserved.
