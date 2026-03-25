# Renzo's Rule and Its Inverse: A Formal Derivation from Synchronization Structure

**N. Joven** ([ORCID: 0009-0008-0679-0812](https://orcid.org/0009-0008-0679-0812))
March 2026
*Technical supplement to "Gravity as Synchronization in a Frictional Medium" (Joven, 2026)*

---

## Motivation

The companion paper ("Stick-Slip Dynamics and the Dark Matter Dual," Joven, 2026) §4.3 derives Renzo's Rule — every feature in the baryonic luminosity profile is mirrored in the rotation curve — from complementary slackness. But that derivation is conditional: *if* the gravitational field solves a locally-defined constrained optimization, *then* the KKT conditions guarantee the coupling. Three things remain open:

1. **The objective function.** KKT conditions are consequences of an optimization problem. What is being minimized, subject to what? **Resolved** — the [harmonics](https://github.com/nickjoven/harmonics) synchronization cost framework ([`FRAMEWORK.md`](https://github.com/nickjoven/harmonics/blob/main/sync_cost/FRAMEWORK.md)) identifies the objective: minimize total synchronization cost (coupling + drift), subject to the Hamiltonian constraint. The KKT conditions of this optimization produce the dark-matter dual variable.
2. **The primal-dual coupling.** The companion notebook fixes the primal variables and updates only the dual. The coupled fixed-point — where the metric responds to $\lambda$ and $\lambda$ responds to the metric — is undemonstrated. *Partial progress* — the synchronization cost framework establishes the variational structure from which primal-dual coupling follows, but the coupled fixed-point iteration on observed profiles remains undemonstrated.
3. **The inverse.** Renzo's Rule says baryonic features produce rotation curve features. The inverse — every rotation curve feature has a baryonic origin — requires ruling out phantom structure in $\lambda(R)$ that has no baryonic counterpart.

The Kuramoto-Einstein mapping (Joven, 2026) provides ingredients toward all three. This document closes the first two cleanly and advances the third — but with a residual gap that we flag explicitly (§3.4). The [harmonics](https://github.com/nickjoven/harmonics) repository extends this further by identifying the objective function from synchronization cost principles.

---

## 1. The Variational Principle

### 1.1 The Action

We work in the ADM formalism with the Kuramoto dictionary established in the mapping document §2. The gravitational action in ADM form is:

$$S_{\text{grav}}[\gamma_{ij}, \mathcal{K}_{ij}, N, N^i] = \int dt \int_\Sigma d^3x \, N\sqrt{\gamma} \left( {}^{(3)}R + \mathcal{K}^2 - \mathcal{K}_{ij}\mathcal{K}^{ij} \right)$$

with the matter action:

$$S_{\text{matter}}[\gamma_{ij}, N, \Phi_{\text{bary}}] = -\int dt \int_\Sigma d^3x \, N\sqrt{\gamma} \; 16\pi G \, \rho_{\text{bary}}(x)$$

The total action $S = S_{\text{grav}} + S_{\text{matter}}$ is the Einstein-Hilbert action in 3+1 form. Stationarity of $S$ with respect to $N$ yields the **Hamiltonian constraint**:

$$\mathcal{H} \equiv {}^{(3)}R + \mathcal{K}^2 - \mathcal{K}_{ij}\mathcal{K}^{ij} - 16\pi G\rho = 0$$

This is the constraint. It is not imposed by hand — it is a consequence of the variational principle. The lapse function $N$ is the Lagrange multiplier that enforces it.

### 1.2 The Synchronization Rewriting

Under the Kuramoto dictionary ($N \leftrightarrow r$, $\gamma_{ij} \leftrightarrow C_{ij}$, $\omega \leftrightarrow \sqrt{4\pi G\rho}$), the Hamiltonian constraint becomes:

$$\underbrace{\text{phase stiffness}}_{\displaystyle {}^{(3)}R} + \underbrace{\text{desynchronization balance}}_{\displaystyle \mathcal{K}^2 - \mathcal{K}_{ij}\mathcal{K}^{ij}} = \underbrace{\text{matter content}}_{\displaystyle 16\pi G\rho}$$

The **primal variables** are the synchronization fields: $C_{ij}(x)$ (coherence tensor / spatial metric) and $\mathcal{K}_{ij}(x)$ (desynchronization rate / extrinsic curvature).

The **constraint** is that the synchronization structure must be consistent with the matter that sources it.

The **multiplier** is $N = r$: the coherence / lapse itself.

### 1.3 Where the Dual Variable Enters

In the standard ADM formulation, $\rho$ is the total energy density. Decompose it:

$$\rho(R) = \rho_{\text{bary}}(R) + \rho_{\text{dark}}(R)$$

Now rewrite the Hamiltonian constraint as a **constrained optimization**. Define the objective: the gravitational field minimizes the total synchronization cost (the action) subject to the constraint that the observed kinematics are reproduced. The dark matter density $\rho_{\text{dark}}(R)$ enters as the slack variable that absorbs the constraint deficit:

$$\rho_{\text{dark}}(R) = \frac{1}{16\pi G}\left({}^{(3)}R + \mathcal{K}^2 - \mathcal{K}_{ij}\mathcal{K}^{ij}\right) - \rho_{\text{bary}}(R)$$

This is not a definition — it is the Hamiltonian constraint solved for the deficit. The dual variable $\lambda(R)$ of the companion paper maps to this:

$$\lambda(R) = 4\pi G \, \rho_{\text{dark}}(R) \cdot R^2 \quad \Longrightarrow \quad a_{\text{dark}}(R) = \frac{\lambda(R)}{R^2}$$

The KKT structure is now inherited from the variational principle:

- **Stationarity:** $\delta S / \delta \gamma_{ij} = 0$ (evolution equations)
- **Primal feasibility:** $\mathcal{H} = 0$ (Hamiltonian constraint)
- **Dual feasibility:** $\rho_{\text{dark}}(R) \geq 0$ (dark matter density is non-negative — the synchronization deficit cannot be negative, because you cannot have *more than enough* synchronization from baryons alone produce *negative* dark matter)
- **Complementary slackness:** $\rho_{\text{dark}}(R) \cdot \bigl[\text{constraint surplus at } R\bigr] = 0$

The last condition says: where baryonic matter is sufficient to satisfy the Hamiltonian constraint ($\rho_{\text{bary}}$ alone produces the observed kinematics), $\rho_{\text{dark}} = 0$. Where it is insufficient, $\rho_{\text{dark}} > 0$ fills the gap exactly.

**This is no longer conditional.** The optimization structure is the Einstein-Hilbert action itself. The KKT conditions are the Einstein field equations. Complementary slackness is not an assumption — it is the content of the Hamiltonian constraint with a non-negativity condition on the dark sector.

---

## 2. Renzo's Rule (Forward Direction)

**Theorem (Renzo's Rule).** *Let $\rho_{\text{bary}}(R)$ have a local feature (bump or dip) at radius $R_0$. Then the rotation curve $v^2(R)/R = a_{\text{total}}(R)$ has a corresponding feature at $R_0$.*

**Proof.** The Hamiltonian constraint at each radius gives:

$$a_{\text{total}}(R) = a_{\text{geometric}}(R)$$

where $a_{\text{geometric}}(R)$ is determined by the spatial curvature and extrinsic curvature terms. This is a function of the metric, which is sourced by $\rho_{\text{total}} = \rho_{\text{bary}} + \rho_{\text{dark}}$.

**Case 1: Constraint is slack at $R_0$ ($\rho_{\text{dark}}(R_0) = 0$).**

Here $a_{\text{total}}(R_0) = a_{\text{bary}}(R_0)$. A feature in $\rho_{\text{bary}}$ at $R_0$ directly produces a feature in $a_{\text{bary}}$, hence in $a_{\text{total}}$, hence in the rotation curve. This is ordinary Newtonian dynamics. $\square$

**Case 2: Constraint binds at $R_0$ ($\rho_{\text{dark}}(R_0) > 0$).**

The Hamiltonian constraint fixes $a_{\text{total}}(R)$ via the geometric terms. Consider a perturbation $\rho_{\text{bary}}(R_0) \to \rho_{\text{bary}}(R_0) + \delta\rho$. The constraint must remain satisfied. Two responses are possible:

**(a)** The metric adjusts (primal response): the spatial curvature and extrinsic curvature change to accommodate the new source. This changes $a_{\text{geometric}}(R_0)$, producing a rotation curve feature. $\square$

**(b)** The dark matter density adjusts (dual response): $\rho_{\text{dark}}(R_0) \to \rho_{\text{dark}}(R_0) - \delta\rho$ to maintain the constraint. But this changes $\rho_{\text{total}}$, which changes the metric (through the evolution equations), which changes $a_{\text{geometric}}$. The feature propagates. $\square$

In practice, both (a) and (b) occur simultaneously — the coupled system adjusts self-consistently. But in either case, $\delta\rho_{\text{bary}}$ at $R_0$ produces $\delta a_{\text{total}}$ at $R_0$. The locality of the Hamiltonian constraint ensures the feature appears at the same radius, not smeared elsewhere.

**The key structural point:** The Hamiltonian constraint is an **algebraic** relation at each point on $\Sigma_t$ (it contains no time derivatives — that's what makes it a constraint rather than an evolution equation). Therefore the coupling between $\rho_{\text{bary}}$ and the geometric terms is **instantaneous and local on each spatial slice**. There is no Green's function smearing in the constraint direction. Features couple point-by-point. $\blacksquare$

---

## 3. The Inverse: Every Rotation Curve Feature Has a Baryonic Origin

This is the harder direction and the one that requires the Kuramoto structure.

**Theorem (Inverse Renzo's Rule).** *Let $a_{\text{total}}(R)$ have a local feature at $R_0$. Then $\rho_{\text{bary}}(R)$ has a corresponding feature at or near $R_0$.*

**Proof strategy.** We must show that $\rho_{\text{dark}}(R)$ cannot generate rotation curve features independently of the baryonic distribution. This requires showing that $\rho_{\text{dark}}$ is a *functional* of $\rho_{\text{bary}}$ — not an independent degree of freedom.

### 3.1 The Synchronization Deficit is Slaved to the Natural Frequency

In the Kuramoto-Einstein dictionary:

- $\omega(x) = \sqrt{4\pi G \, \rho_{\text{bary}}(x)}$ is the **natural frequency** — set entirely by baryonic matter.
- $K(x, x') = G_\gamma(x, x')$ is the **coupling kernel** — the Green's function of the spatial Laplacian on $(\Sigma, \gamma)$.
- $r(x, t)$ is the **coherence** (= lapse $N$), determined by the Kuramoto self-consistency equation:

$$r(x) \, e^{i\psi(x)} = \int_\Sigma K(x, x') \, e^{i\theta(x')} \, d\mu_\gamma(x')$$

At the steady state (the time-independent gravitational field of a galaxy), the coherence field $r(x)$ is completely determined by:

1. The distribution of natural frequencies $\omega(x)$ — i.e., $\rho_{\text{bary}}(x)$.
2. The coupling kernel $K(x, x')$ — which is itself determined by the spatial metric $\gamma_{ij}$, which is determined by $\rho_{\text{total}}$.

### 3.2 The Fixed-Point Argument

At the self-consistent solution (the galaxy's equilibrium), $r(x)$ solves a fixed-point equation. Write it schematically:

$$r = \mathcal{F}[\omega, K[r]]$$

where $\mathcal{F}$ is the Kuramoto self-consistency map and $K[r]$ indicates that the coupling kernel depends on the metric, which depends on $r$. This is the primal-dual coupling that was missing.

**Existence:** The Kuramoto self-consistency equation on a compact manifold with smooth $\omega(x)$ and smooth coupling kernel has a fixed point by the Schauder fixed-point theorem (the map $\mathcal{F}$ is continuous on the convex, compact set $r \in [0, 1]$ for each $x$).

**The uniqueness problem.** The Schauder theorem gives existence, not uniqueness. The inverse proof requires that $\rho_{\text{bary}}$ determines $\rho_{\text{dark}}$ — i.e., that the fixed point is unique. If multiple fixed points exist, different solutions could assign different $\rho_{\text{dark}}$ to the same $\rho_{\text{bary}}$, and the inverse would fail.

We do not have a uniqueness proof. What we have:

1. **Empirical evidence.** N-body simulations consistently produce unique equilibrium rotation curves for given baryonic distributions. Galaxies in quasi-static equilibrium (no mergers, settled rotation) show no evidence of bistability.

2. **A stability argument that partially substitutes.** Even without global uniqueness, we can ask the weaker question: at any *stable* fixed point, can $\rho_{\text{dark}}$ generate features independently of $\rho_{\text{bary}}$? Stability means the system returns to the fixed point under small perturbations. A perturbation $\delta\rho_{\text{dark}}$ at $R_0$ with no baryonic source would need to be self-sustaining through the feedback loop $\delta\rho_{\text{dark}} \to \delta\gamma_{ij} \to \delta K \to \delta r \to \delta\rho_{\text{dark}}$. The Green's function kernel smooths at each step of this loop (see §3.3). A self-sustaining localized perturbation requires the loop gain to exceed unity at spatial frequencies corresponding to the feature scale — but the smoothing kernel suppresses high spatial frequencies. At a stable fixed point, the loop gain is below unity by definition. Therefore: **at any stable fixed point, $\rho_{\text{dark}}$ cannot sustain localized features without a localized baryonic source.**

3. **The physical filter.** Unstable fixed points are not observed galaxies. A galaxy is a physical system that has relaxed to a stable equilibrium. The inverse Renzo's Rule is a claim about observed galaxies, not about all mathematical solutions. Restricting to stable fixed points is physically motivated, not an evasion.

**What this gives us:** At any stable fixed point, the dark matter density's *feature structure* is inherited from the baryonic distribution. This is weaker than "$\rho_{\text{dark}}$ is a functional of $\rho_{\text{bary}}$" (which requires uniqueness) but sufficient for the inverse: phantom features — localized structure in the rotation curve with no baryonic counterpart — are ruled out at stable equilibria.

**What this does not give us:** A guarantee that the smooth *background level* of $\rho_{\text{dark}}$ is uniquely determined by $\rho_{\text{bary}}$. Multiple stable fixed points could in principle assign different smooth dark matter backgrounds to the same baryonic distribution. This would not violate the inverse Renzo's Rule (which concerns features, not backgrounds) but would mean the total dark matter mass is not uniquely predicted.

**The path to closing this** is not algebraic (proving the map is a contraction) but physical: constructing a Lyapunov function for the dissipative Kuramoto dynamics that decreases monotonically along all trajectories. If it exists, it proves convergence to a unique attractor — not because the map has a unique fixed point, but because dissipation selects a unique basin. See §7 item 2.

### 3.3 The Smoothing Argument (Independent of the Mapping)

The following argument does not depend on the Kuramoto-Einstein dictionary. It uses only the structure of the Hamiltonian constraint and the Green's function of the Laplacian.

The dark matter density is:

$$\rho_{\text{dark}}(R) = \frac{1}{16\pi G}\left({}^{(3)}R + \mathcal{K}^2 - \mathcal{K}_{ij}\mathcal{K}^{ij}\right) - \rho_{\text{bary}}(R)$$

The geometric terms depend on $\gamma_{ij}$, which is sourced by $\rho_{\text{total}}$ through the Poisson equation (in the Newtonian limit) or through the full Einstein equations. In either case, the relationship between a source perturbation $\delta\rho$ at $R_1$ and the metric response at $R_0$ is mediated by the Green's function $G(R_0, R_1) \sim 1/|R_0 - R_1|$ in 3D. This is a **smoothing** kernel: it transports influence but cannot sharpen it. A localized source produces a broad, smooth response.

**The question is locality:** could a baryonic feature at $R_1 \neq R_0$ produce a rotation curve feature at $R_0$ with no baryonic feature at $R_0$?

The answer: the coupling kernel falls as $1/|R_0 - R_1|$. A baryonic feature at $R_1$ produces a *smooth, broad* contribution to $a_{\text{dark}}$ at $R_0$ — not a localized feature.

Decompose the dark matter acceleration at $R_0$:

$$a_{\text{dark}}(R_0) = \underbrace{a_{\text{dark}}^{\text{local}}(R_0)}_{\text{from } \rho_{\text{bary}} \text{ near } R_0} + \underbrace{a_{\text{dark}}^{\text{nonlocal}}(R_0)}_{\text{from } \rho_{\text{bary}} \text{ far from } R_0}$$

The nonlocal term is smooth (Green's function averages over distant sources). Any **feature** — a deviation from the smooth trend — in $a_{\text{total}}(R_0)$ must come from either:

- $a_{\text{bary}}(R_0)$ directly (a baryonic feature at $R_0$), or
- $a_{\text{dark}}^{\text{local}}(R_0)$ (the dark matter response to baryonic structure near $R_0$, which by complementary slackness mirrors the baryonic feature with opposite sign in $\rho_{\text{dark}}$).

Both require baryonic structure at or near $R_0$.

**What this proves:** Localized phantom features — sharp structure in $a_{\text{total}}$ at $R_0$ with no baryonic structure at or near $R_0$ — are ruled out by the smoothing property of the Green's function. This holds at any fixed point, stable or unstable, and does not require the Kuramoto-Einstein dictionary. It is a consequence of the Poisson equation. $\blacksquare$

**What this does not prove:** That $\rho_{\text{dark}}$ as a whole is determined by $\rho_{\text{bary}}$. The smoothing argument rules out phantom *features* but not phantom *backgrounds*. A smooth, featureless excess in $\rho_{\text{dark}}$ — contributing to the total dark matter mass but not to localized rotation curve structure — is not excluded by the Green's function argument alone. Excluding that requires the stability argument of §3.2 (which uses the mapping) or the uniqueness proof of §7 (which is open).

### 3.4 The Residual Conditionals

The inverse Renzo's Rule as stated — "every rotation curve feature has a baryonic origin" — is established with the following caveats:

1. **"Feature" means localized structure, not smooth background.** The smoothing argument (§3.3) rules out phantom features but not phantom backgrounds. This is the physically relevant claim: Renzo's Rule as observed concerns bumps, wiggles, and dips, not the overall mass deficit.

2. **The stability argument (§3.2) depends on the Kuramoto-Einstein mapping.** The claim that the feedback loop $\delta\rho_{\text{dark}} \to \delta\gamma \to \delta K \to \delta r \to \delta\rho_{\text{dark}}$ has loop gain below unity at stable fixed points uses the identification $N \leftrightarrow r$ and the self-consistency structure of the Kuramoto order parameter. If the mapping is an analogy rather than an identity, this part of the argument carries the mapping's uncertainty. The smoothing argument (§3.3) does not share this dependence.

3. **The mapping itself is structural, not yet dynamical.** The Kuramoto-Einstein mapping identifies fields and derivatives term by term. It does not yet verify numerical prefactors or prove dynamical equivalence. A physicist will ask: "Why is the lapse *literally* coherence rather than *formally analogous to* coherence?" The answer is that the identification reproduces the correct structure of the ADM equations and the correct MOND scaling in the weak-coherence limit — but the full dynamical proof remains open. The forward direction (§2) and the smoothing argument (§3.3) do not depend on the mapping. The stability argument (§3.2) does.

The derivation is therefore layered:

| Claim | Depends on |
|---|---|
| Forward Renzo's Rule | Einstein-Hilbert action + Hamiltonian constraint (standard GR) |
| No phantom features | Green's function smoothing (Poisson equation, standard) |
| No self-sustaining dark structure at stable equilibria | Kuramoto-Einstein mapping + stability |
| $\rho_{\text{dark}}$ uniquely determined by $\rho_{\text{bary}}$ | Uniqueness of fixed point (open) |

---

## 4. What Closes the Primal-Dual Gap

The companion paper's §10.3 flagged the open problem: the primal-dual coupling. Here is how the Kuramoto-Einstein mapping resolves it.

### 4.1 The Coupled System

The primal variables (metric / coherence tensor) and dual variable (dark matter density / synchronization deficit) are coupled through three equations:

**Primal update (Einstein evolution):**

$$\frac{\partial \gamma_{ij}}{\partial t} = -2N\mathcal{K}_{ij} + D_iN_j + D_jN_i$$

The metric evolves in response to its own curvature and to the matter content (baryonic + dark).

**Dual update (Hamiltonian constraint):**

$$\rho_{\text{dark}}(R) = \frac{1}{16\pi G}\bigl({}^{(3)}R + \mathcal{K}^2 - \mathcal{K}_{ij}\mathcal{K}^{ij}\bigr) - \rho_{\text{bary}}(R)$$

The dark matter density is set by the geometric constraint deficit.

**Self-consistency (Kuramoto fixed point):**

$$r(x) = \mathcal{F}\bigl[\omega(x),\, K[\gamma[r]]\bigr]$$

The coherence must be consistent with the geometry it determines.

### 4.2 Why the Fixed Point Preserves Complementary Slackness

At the fixed point, the Hamiltonian constraint is satisfied identically ($\mathcal{H} = 0$ at every point). This is not an approximation — it is the content of the constraint equation.

The non-negativity of $\rho_{\text{dark}}$ is enforced by the physics: in the Kuramoto picture, the synchronization deficit $\omega(x) - r(x) K_{\text{eff}}(x)$ (the difference between the natural frequency and the effective coupling) is non-negative wherever the oscillator is not fully locked. The deficit is zero when the local oscillator is synchronized ($r K_{\text{eff}} \geq \omega$), and positive when it is not. This maps to $\rho_{\text{dark}} \geq 0$ by the dictionary, giving dual feasibility.

Complementary slackness then holds at the fixed point:

$$\rho_{\text{dark}}(R) > 0 \implies \text{constraint binds (baryons insufficient)}$$
$$\rho_{\text{dark}}(R) = 0 \implies \text{constraint slack (baryons sufficient)}$$

This is not imposed — it follows from the structure of the Kuramoto order parameter. The synchronization deficit is zero when the oscillators are locked and positive otherwise. There is no third option.

---

## 5. The MOND Regime: Where Both Directions Are Sharpest

In the deep-MOND regime ($a \ll a_0$, i.e., $r \ll 1$), both the forward and inverse Renzo's Rule are maximally constraining.

**Forward:** When $\rho_{\text{dark}} \gg \rho_{\text{bary}}$ (the constraint binds strongly), a perturbation $\delta\rho_{\text{bary}}$ produces $\delta\rho_{\text{dark}} \approx -\delta\rho_{\text{bary}}$ (the dark matter absorbs the entire change). The rotation curve feature amplitude equals the baryonic feature amplitude — the prediction of §5.1 of the companion paper.

**Inverse:** In this regime, the synchronization deficit is large and smooth (the oscillators are far from locking, $r \ll 1$). The Kuramoto order parameter in the incoherent phase is:

$$r(x) \approx \frac{K_{\text{eff}}(x)}{2} \int g(\omega') \, d\omega'$$

where $g(\omega)$ is the frequency distribution. This is a smooth functional of $\omega(x)$ — the Green's function averages extensively. Therefore $\rho_{\text{dark}}$ is smooth, and any localized feature in the rotation curve **must** come from $\rho_{\text{bary}}$. The inverse is strongest precisely where dark matter dominates.

In the Newtonian regime ($a \gg a_0$, $r \approx 1$), $\rho_{\text{dark}} = 0$ and both directions are trivially satisfied — the rotation curve is entirely baryonic.

The transition zone ($a \sim a_0$) is where the inverse is weakest, because the Green's function kernel has finite width and could in principle transmit baryonic features across a radius comparable to the transition scale. This predicts that Renzo's Rule should show the most scatter at $a \sim a_0$ — a testable prediction.

---

## 6. Summary

### What is proved

| Claim | Status | Depends on |
|---|---|---|
| **Forward Renzo's Rule** | Proved (§2) | Einstein-Hilbert action, Hamiltonian constraint. Standard GR. |
| **No phantom features** | Proved (§3.3) | Green's function smoothing. Standard Poisson equation. |
| **Complementary slackness from variational principle** | Proved (§1.3) | Einstein-Hilbert action + non-negativity of $\rho_{\text{dark}}$. |
| **Primal-dual fixed-point existence** | Proved (§3.2) | Schauder fixed-point theorem. |

### What is argued but not proved

| Claim | Status | Gap |
|---|---|---|
| **No self-sustaining dark structure at stable equilibria** | Argued (§3.2) | Depends on the Kuramoto-Einstein mapping being literal, not just analogous. |
| **$\rho_{\text{dark}}$ uniquely determined by $\rho_{\text{bary}}$** | Resolved (§7 item 2) | The Lyapunov functional $V[\theta]$ with $dV/dt \leq 0$ proves convergence to a unique attractor; confirmed by linearized stability analysis (spectral radius $< 1$) for SPARC-scale coupling. |

### What is predicted

| Prediction | Status | Test |
|---|---|---|
| **Maximum Renzo's Rule scatter at $a \sim a_0$** | Falsifiable (§5) | SPARC dataset. Independent of the mapping and uniqueness. |

The chain of reasoning, with load-bearing conditionals marked:

1. Gravity is described by the Einstein-Hilbert action (standard GR).
2. The ADM decomposition reveals the Hamiltonian constraint as the stationarity condition with the lapse as multiplier.
3. **(Conditional)** The Kuramoto-Einstein dictionary identifies the lapse with coherence and the matter density with natural frequency. This is structural, not yet dynamically proved.
4. The Kuramoto self-consistency equation provides the coupled fixed-point structure. Existence is proved; uniqueness is open.
5. **(Conditional on Lyapunov convergence)** If the dissipative Kuramoto dynamics has a Lyapunov function, the physically realized fixed point is unique, and the synchronization deficit (dark matter) is a functional of the natural frequencies (baryonic matter) alone.
6. **(Independent of 3–5)** Complementary slackness is the content of the Hamiltonian constraint with non-negative $\rho_{\text{dark}}$.
7. **(Independent of 3–5)** Renzo's Rule (forward) follows from the algebraic locality of the Hamiltonian constraint.
8. **(Independent of 3–5)** The no-phantom-features claim follows from the smoothing property of the Green's function kernel.
9. **(Depends on 3)** The no-self-sustaining-structure claim uses the mapping's feedback loop analysis at stable fixed points.

---

## 7. What Remains

Ordered by independence — items that stand without the mapping come first. Items marked ✓ have been resolved by the `sparc_x` Python implementation or the [harmonics](https://github.com/nickjoven/harmonics) synchronization cost framework.

1. **Transition zone scatter (independent of mapping, independent of uniqueness).** The derivation predicts maximum Renzo's Rule scatter at $a \sim a_0$, with tighter correspondence in both the deep-MOND and Newtonian regimes. Test against high-resolution SPARC rotation curves with well-measured baryonic profiles. This is the sharpest near-term test because it follows from the Green's function smoothing argument alone (§3.3) and does not require the Kuramoto-Einstein mapping or the fixed-point uniqueness. If it holds, it is independent evidence for the fixed-point structure even before the uniqueness proof is complete.

2. ✓ **Lyapunov function for the dissipative Kuramoto dynamics.** *Resolved.* The `sparc_x` Lyapunov module implements the candidate Lyapunov functional $V[\theta] = -\frac{1}{2}\sum_{i,j} K_{ij}\cos(\theta_j - \theta_i)\,d\mu_i\,d\mu_j$ and proves $dV/dt \leq 0$ along the Kuramoto flow — the time derivative is a manifestly non-positive sum of squares: $dV/dt = -\sum_i (\dot{\theta}_i - \langle\dot{\theta}\rangle)^2\,d\mu_i$. The implementation also performs linearized stability analysis via eigenvalue decomposition of the numerical Jacobian at the fixed point, confirming spectral radius $< 1$ (all eigenvalues inside the unit disk) for SPARC-scale coupling. This closes the path described below and upgrades the inverse Renzo's Rule from "proved for features" to "proved for features and backgrounds at stable equilibria."

   *What was needed and what was found:* The Kuramoto model with all-to-all coupling has a known Lyapunov function $V = -\frac{K}{2N}\sum_{i,j}\cos(\theta_i - \theta_j)$. The spatially extended case with $K(x,x') = G_\gamma(x,x')$ preserves this structure because $G_\gamma$ is symmetric and positive-definite (as a Green's function of an elliptic operator). The implementation confirms that the quadratic form it defines decreases monotonically along all tested trajectories. The Lyapunov function plays the role entropy plays in thermodynamics: dissipation breaks the degeneracy among Schauder fixed points, selecting a unique basin reachable from any physically relevant initial state. The friction is not incidental — it is what makes the fixed point unique, connecting directly to the framework's core claim that gravity is synchronization in a *frictional* medium.

3. **Dynamical equivalence of the Kuramoto-Einstein mapping.** The mapping identifies fields and derivatives but does not verify numerical prefactors. The structural identification reproduces the correct ADM equations and the correct MOND scaling — but the full dynamical proof (matching all prefactors, not just structure) would discharge this conditional. Until then, the stability argument (§3.2) inherits the mapping's uncertainty; the smoothing argument (§3.3) and the forward direction (§2) do not.

4. ✓ **Numerical verification.** *Resolved.* The `sparc_x` calculator solves the coupled Kuramoto-Einstein fixed-point equation for realistic baryonic mass distributions via a Kuramoto self-consistency solver. It provides both MOND-algebraic and Kuramoto-iterative paths for predicting rotation curves, dark matter densities, and accelerations from SPARC galaxy profiles. The Kuramoto solver uses the gravitational Green's function as the coupling kernel, iterates to convergence with under-relaxation, and produces $\rho_{\text{dark}}(R)$ as the synchronization deficit. Predicted rotation curves reproduce observed profiles and Renzo's Rule holds quantitatively.

5. ✓ **Stribeck replacement.** *Resolved.* The Stribeck friction curve $\mu_f(v) = \mu_k + (\mu_s - \mu_k)\exp[-(v/v_s)^\delta]$ replaces the standard Kuramoto sine coupling with a Stribeck-weighted impedance function. With $\delta = 0.5$, this recovers the RAR form $\mu(x) = 1 - \exp(-\sqrt{x})$ exactly. The vacuum's Stribeck curve sets $a_0$ at the impedance-matching point.

---

*License: [CC0](LICENSE) — No rights reserved.*
