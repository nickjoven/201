# Renzo's Rule and Its Inverse: A Formal Derivation from Synchronization Structure

**N. Joven** ([ORCID: 0009-0008-0679-0812](https://orcid.org/0009-0008-0679-0812))
March 2026
*Technical supplement to [Gravity as Synchronization in a Frictional Medium](joven_unifying_framework.md)*

---

## Motivation

The [companion paper](intersections/joven_stick_slip_dark_matter.md) §4.3 derives Renzo's Rule — every feature in the baryonic luminosity profile is mirrored in the rotation curve — from complementary slackness. But that derivation is conditional: *if* the gravitational field solves a locally-defined constrained optimization, *then* the KKT conditions guarantee the coupling. Three things remain open:

1. **The objective function.** KKT conditions are consequences of an optimization problem. What is being minimized, subject to what?
2. **The primal-dual coupling.** The companion notebook fixes the primal variables and updates only the dual. The coupled fixed-point — where the metric responds to $\lambda$ and $\lambda$ responds to the metric — is undemonstrated.
3. **The inverse.** Renzo's Rule says baryonic features produce rotation curve features. The inverse — every rotation curve feature has a baryonic origin — requires ruling out phantom structure in $\lambda(r)$ that has no baryonic counterpart.

The [Kuramoto-Einstein mapping](kuramoto_einstein_mapping.md) provides the missing ingredients. This document closes all three gaps.

---

## 1. The Variational Principle

### 1.1 The Action

We work in the ADM formalism with the Kuramoto dictionary established in the [mapping document](kuramoto_einstein_mapping.md) §2. The gravitational action in ADM form is:

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

$$\rho(r) = \rho_{\text{bary}}(r) + \rho_{\text{dark}}(r)$$

Now rewrite the Hamiltonian constraint as a **constrained optimization**. Define the objective: the gravitational field minimizes the total synchronization cost (the action) subject to the constraint that the observed kinematics are reproduced. The dark matter density $\rho_{\text{dark}}(r)$ enters as the slack variable that absorbs the constraint deficit:

$$\rho_{\text{dark}}(r) = \frac{1}{16\pi G}\left({}^{(3)}R + \mathcal{K}^2 - \mathcal{K}_{ij}\mathcal{K}^{ij}\right) - \rho_{\text{bary}}(r)$$

This is not a definition — it is the Hamiltonian constraint solved for the deficit. The dual variable $\lambda(r)$ of the companion paper maps to this:

$$\lambda(r) = 4\pi G \, \rho_{\text{dark}}(r) \cdot r^2 \quad \Longrightarrow \quad a_{\text{dark}}(r) = \frac{\lambda(r)}{r^2}$$

The KKT structure is now inherited from the variational principle:

- **Stationarity:** $\delta S / \delta \gamma_{ij} = 0$ (evolution equations)
- **Primal feasibility:** $\mathcal{H} = 0$ (Hamiltonian constraint)
- **Dual feasibility:** $\rho_{\text{dark}}(r) \geq 0$ (dark matter density is non-negative — the synchronization deficit cannot be negative, because you cannot have *more than enough* synchronization from baryons alone produce *negative* dark matter)
- **Complementary slackness:** $\rho_{\text{dark}}(r) \cdot \bigl[\text{constraint surplus at } r\bigr] = 0$

The last condition says: where baryonic matter is sufficient to satisfy the Hamiltonian constraint ($\rho_{\text{bary}}$ alone produces the observed kinematics), $\rho_{\text{dark}} = 0$. Where it is insufficient, $\rho_{\text{dark}} > 0$ fills the gap exactly.

**This is no longer conditional.** The optimization structure is the Einstein-Hilbert action itself. The KKT conditions are the Einstein field equations. Complementary slackness is not an assumption — it is the content of the Hamiltonian constraint with a non-negativity condition on the dark sector.

---

## 2. Renzo's Rule (Forward Direction)

**Theorem (Renzo's Rule).** *Let $\rho_{\text{bary}}(r)$ have a local feature (bump or dip) at radius $r_0$. Then the rotation curve $v^2(r)/r = a_{\text{total}}(r)$ has a corresponding feature at $r_0$.*

**Proof.** The Hamiltonian constraint at each radius gives:

$$a_{\text{total}}(r) = a_{\text{geometric}}(r)$$

where $a_{\text{geometric}}(r)$ is determined by the spatial curvature and extrinsic curvature terms. This is a function of the metric, which is sourced by $\rho_{\text{total}} = \rho_{\text{bary}} + \rho_{\text{dark}}$.

**Case 1: Constraint is slack at $r_0$ ($\rho_{\text{dark}}(r_0) = 0$).**

Here $a_{\text{total}}(r_0) = a_{\text{bary}}(r_0)$. A feature in $\rho_{\text{bary}}$ at $r_0$ directly produces a feature in $a_{\text{bary}}$, hence in $a_{\text{total}}$, hence in the rotation curve. This is ordinary Newtonian dynamics. $\square$

**Case 2: Constraint binds at $r_0$ ($\rho_{\text{dark}}(r_0) > 0$).**

The Hamiltonian constraint fixes $a_{\text{total}}(r)$ via the geometric terms. Consider a perturbation $\rho_{\text{bary}}(r_0) \to \rho_{\text{bary}}(r_0) + \delta\rho$. The constraint must remain satisfied. Two responses are possible:

**(a)** The metric adjusts (primal response): the spatial curvature and extrinsic curvature change to accommodate the new source. This changes $a_{\text{geometric}}(r_0)$, producing a rotation curve feature. $\square$

**(b)** The dark matter density adjusts (dual response): $\rho_{\text{dark}}(r_0) \to \rho_{\text{dark}}(r_0) - \delta\rho$ to maintain the constraint. But this changes $\rho_{\text{total}}$, which changes the metric (through the evolution equations), which changes $a_{\text{geometric}}$. The feature propagates. $\square$

In practice, both (a) and (b) occur simultaneously — the coupled system adjusts self-consistently. But in either case, $\delta\rho_{\text{bary}}$ at $r_0$ produces $\delta a_{\text{total}}$ at $r_0$. The locality of the Hamiltonian constraint ensures the feature appears at the same radius, not smeared elsewhere.

**The key structural point:** The Hamiltonian constraint is an **algebraic** relation at each point on $\Sigma_t$ (it contains no time derivatives — that's what makes it a constraint rather than an evolution equation). Therefore the coupling between $\rho_{\text{bary}}$ and the geometric terms is **instantaneous and local on each spatial slice**. There is no Green's function smearing in the constraint direction. Features couple point-by-point. $\blacksquare$

---

## 3. The Inverse: Every Rotation Curve Feature Has a Baryonic Origin

This is the harder direction and the one that requires the Kuramoto structure.

**Theorem (Inverse Renzo's Rule).** *Let $a_{\text{total}}(r)$ have a local feature at $r_0$. Then $\rho_{\text{bary}}(r)$ has a corresponding feature at or near $r_0$.*

**Proof strategy.** We must show that $\rho_{\text{dark}}(r)$ cannot generate rotation curve features independently of the baryonic distribution. This requires showing that $\rho_{\text{dark}}$ is a *functional* of $\rho_{\text{bary}}$ — not an independent degree of freedom.

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

**Uniqueness in the galactic regime:** For galaxies in quasi-static equilibrium (no mergers, settled rotation), the system is well past the transient. Numerical N-body simulations consistently produce unique equilibrium rotation curves for given baryonic distributions, providing empirical evidence that the fixed point is unique in this regime.

**The critical consequence:** At the fixed point, $r(x)$ is a functional of $\omega(x)$ alone:

$$r(x) = r[\omega](x) = r[\rho_{\text{bary}}](x)$$

Because the coupling kernel at equilibrium is itself determined by the solution, and the solution is unique, the entire gravitational field — including the dark matter density — is determined by the baryonic distribution.

### 3.3 Ruling Out Phantom Features

Now we can prove the inverse. The dark matter density is:

$$\rho_{\text{dark}}(r) = \frac{1}{16\pi G}\left({}^{(3)}R + \mathcal{K}^2 - \mathcal{K}_{ij}\mathcal{K}^{ij}\right) - \rho_{\text{bary}}(r)$$

The geometric terms on the right depend on $\gamma_{ij}$, which depends on $\rho_{\text{total}}$, which at the fixed point is a functional of $\rho_{\text{bary}}$ alone (§3.2). Therefore:

$$\rho_{\text{dark}}(r) = \rho_{\text{dark}}[\rho_{\text{bary}}](r)$$

$\rho_{\text{dark}}$ is not an independent field. It is **determined** by $\rho_{\text{bary}}$ through the fixed-point structure of the Kuramoto-Einstein equations.

**Now consider a feature in $a_{\text{total}}$ at $r_0$.** Since $a_{\text{total}} = a_{\text{bary}} + a_{\text{dark}}$ and both terms are functionals of $\rho_{\text{bary}}$, the feature in $a_{\text{total}}$ must originate from $\rho_{\text{bary}}$.

**The question is locality:** could a baryonic feature at $r_1 \neq r_0$ produce a rotation curve feature at $r_0$ with no baryonic feature at $r_0$? This would violate the inverse in spirit, though not in letter.

The answer is: the coupling kernel $K(x, x') = G_\gamma(x, x')$ is the Green's function of the Laplacian, which in 3D falls as $1/|x - x'|$. It is a **smoothing** kernel. It can transport influence from $r_1$ to $r_0$, but it does so by spreading, not by sharpening. A baryonic feature at $r_1$ produces a *smooth, broad* contribution to $a_{\text{dark}}$ at $r_0$ — not a localized feature.

Therefore: **a localized feature in $a_{\text{total}}$ at $r_0$ requires a localized source at $r_0$.** The only localized source is $\rho_{\text{bary}}(r_0)$, because the dark matter contribution at $r_0$ from distant baryonic features is smoothed by the Green's function.

More precisely, decompose the dark matter acceleration at $r_0$:

$$a_{\text{dark}}(r_0) = \underbrace{a_{\text{dark}}^{\text{local}}(r_0)}_{\text{from } \rho_{\text{bary}} \text{ near } r_0} + \underbrace{a_{\text{dark}}^{\text{nonlocal}}(r_0)}_{\text{from } \rho_{\text{bary}} \text{ far from } r_0}$$

The nonlocal term is smooth (Green's function averages over distant sources). Any **feature** — a deviation from the smooth trend — in $a_{\text{total}}(r_0)$ must come from either:

- $a_{\text{bary}}(r_0)$ directly (a baryonic feature at $r_0$), or
- $a_{\text{dark}}^{\text{local}}(r_0)$ (the dark matter response to baryonic structure near $r_0$, which by complementary slackness mirrors the baryonic feature with opposite sign in $\rho_{\text{dark}}$).

Both require baryonic structure at or near $r_0$. $\blacksquare$

---

## 4. What Closes the Primal-Dual Gap

The companion paper's §10.3 flagged the open problem: the primal-dual coupling. Here is how the Kuramoto-Einstein mapping resolves it.

### 4.1 The Coupled System

The primal variables (metric / coherence tensor) and dual variable (dark matter density / synchronization deficit) are coupled through three equations:

**Primal update (Einstein evolution):**

$$\frac{\partial \gamma_{ij}}{\partial t} = -2N\mathcal{K}_{ij} + D_iN_j + D_jN_i$$

The metric evolves in response to its own curvature and to the matter content (baryonic + dark).

**Dual update (Hamiltonian constraint):**

$$\rho_{\text{dark}}(r) = \frac{1}{16\pi G}\bigl({}^{(3)}R + \mathcal{K}^2 - \mathcal{K}_{ij}\mathcal{K}^{ij}\bigr) - \rho_{\text{bary}}(r)$$

The dark matter density is set by the geometric constraint deficit.

**Self-consistency (Kuramoto fixed point):**

$$r(x) = \mathcal{F}\bigl[\omega(x),\, K[\gamma[r]]\bigr]$$

The coherence must be consistent with the geometry it determines.

### 4.2 Why the Fixed Point Preserves Complementary Slackness

At the fixed point, the Hamiltonian constraint is satisfied identically ($\mathcal{H} = 0$ at every point). This is not an approximation — it is the content of the constraint equation.

The non-negativity of $\rho_{\text{dark}}$ is enforced by the physics: in the Kuramoto picture, the synchronization deficit $\omega(x) - r(x) K_{\text{eff}}(x)$ (the difference between the natural frequency and the effective coupling) is non-negative wherever the oscillator is not fully locked. The deficit is zero when the local oscillator is synchronized ($r K_{\text{eff}} \geq \omega$), and positive when it is not. This maps to $\rho_{\text{dark}} \geq 0$ by the dictionary, giving dual feasibility.

Complementary slackness then holds at the fixed point:

$$\rho_{\text{dark}}(r) > 0 \implies \text{constraint binds (baryons insufficient)}$$
$$\rho_{\text{dark}}(r) = 0 \implies \text{constraint slack (baryons sufficient)}$$

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

| Gap | Resolution |
|---|---|
| **Objective function unspecified** | The objective is the Einstein-Hilbert action in ADM form. The Hamiltonian constraint is the KKT stationarity condition with $N$ as the multiplier. |
| **Primal-dual coupling open** | The Kuramoto self-consistency equation provides the fixed-point structure. Existence follows from Schauder; uniqueness in the galactic regime from empirical convergence of N-body simulations. |
| **Inverse unproved** | At the fixed point, $\rho_{\text{dark}}$ is a functional of $\rho_{\text{bary}}$. The Green's function kernel is smoothing, so localized rotation curve features require localized baryonic sources. |
| **Complementary slackness assumed** | Complementary slackness is the Kuramoto order parameter's zero/nonzero structure mapped through the dictionary. It is a property of synchronization, not an assumption. |

The chain of reasoning:

1. Gravity is described by the Einstein-Hilbert action (standard GR).
2. The ADM decomposition reveals the Hamiltonian constraint as the stationarity condition with the lapse as multiplier.
3. The Kuramoto-Einstein dictionary identifies the lapse with coherence and the matter density with natural frequency.
4. The Kuramoto self-consistency equation provides the coupled fixed-point structure.
5. At the fixed point, the synchronization deficit (dark matter) is a functional of the natural frequencies (baryonic matter) alone.
6. Complementary slackness is the order parameter's threshold behavior: zero below synchronization, positive above.
7. Renzo's Rule (forward) follows from the algebraic locality of the Hamiltonian constraint.
8. The inverse follows from the smoothing property of the Green's function kernel: $\rho_{\text{dark}}$ inherits its structure from $\rho_{\text{bary}}$ and cannot generate localized features independently.

---

## 7. What Remains

1. **Numerical verification.** Solve the coupled Kuramoto-Einstein fixed-point equation for a realistic baryonic mass distribution (e.g., SPARC galaxies) and verify that the resulting $\rho_{\text{dark}}(r)$ reproduces observed rotation curves with Renzo's Rule holding quantitatively.

2. **Uniqueness proof.** Replace the empirical argument for fixed-point uniqueness with a mathematical proof. Candidate approach: show that the Kuramoto self-consistency map is a contraction in an appropriate function space when $\omega(x)$ has the radial monotonicity typical of galaxies.

3. **Transition zone scatter.** The derivation predicts maximum Renzo's Rule scatter at $a \sim a_0$. Test against high-resolution SPARC rotation curves with well-measured baryonic profiles.

4. **Stribeck replacement.** The standard Kuramoto coupling uses $\sin(\psi - \theta)$. The framework claims the vacuum has a Stribeck-type impedance. Replacing the sine coupling with a Stribeck-weighted function and re-deriving the fixed-point structure would tighten the connection to the Universal Rosin.

---

*License: [CC0](LICENSE) — No rights reserved.*
