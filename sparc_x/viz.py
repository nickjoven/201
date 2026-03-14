"""
Interactive visualizations for the synchronization-gravity framework.

Renders the Kuramoto computation DAG as linked Plotly traces:
each node in the pipeline (baryonic profile -> a_N -> x -> mu -> a_obs ->
v_pred -> deficit -> rho_dark) is an inspectable, hoverable trace.

Requires: plotly (optional dependency).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from sparc_x.calculator import Calculator
from sparc_x.constants import A0, KPC, KMS
from sparc_x.kuramoto import KuramotoState
from sparc_x.stribeck import stribeck_curve, mond_from_stribeck


# ---------------------------------------------------------------------------
# Serialization helpers
# ---------------------------------------------------------------------------

def calculator_to_dict(calc: Calculator) -> dict[str, Any]:
    """Serialize all computed quantities to a JSON-friendly dict.

    This is the bridge between the Python computation and any frontend.
    Every array becomes a list; every float stays a float.
    """
    p = calc.profile
    masks = calc.get_regime_mask()
    accel = calc.get_accelerations()

    d: dict[str, Any] = {
        "name": p.name,
        "r_kpc": p.r_kpc.tolist(),
        "v_obs": p.v_obs.tolist(),
        "e_v_obs": p.e_v_obs.tolist(),
        "v_gas": p.v_gas.tolist(),
        "v_disk": p.v_disk.tolist(),
        "v_bul": p.v_bul.tolist(),
        # Kuramoto DAG nodes
        "a_N": accel["a_N"].tolist(),
        "a_obs": accel["a_obs"].tolist(),
        "a_pred": accel["a_pred"].tolist(),
        "x": p.x.tolist(),
        "v_mond": calc.get_rotation_curve(method="mond").tolist(),
        "rho_dark": calc.get_dark_matter_density().tolist(),
        "prefactors": calc.get_prefactors().tolist(),
        # regime masks
        "regime": [
            "newtonian" if masks["newtonian"][i]
            else ("deep_mond" if masks["deep_mond"][i] else "transition")
            for i in range(len(p.r_kpc))
        ],
        "a0": calc.a0,
        "interpolation": calc.interpolation,
    }

    # Kuramoto state if solved
    if calc._results.kuramoto_state is not None:
        ks = calc._results.kuramoto_state
        d["kuramoto"] = {
            "coherence": ks.coherence.tolist(),
            "mean_phase": ks.mean_phase.tolist(),
            "K_eff": ks.K_eff.tolist(),
            "n_iter": ks.n_iter,
        }
        d["v_kuramoto"] = calc.get_rotation_curve(method="kuramoto").tolist()

    # Lyapunov if computed
    if calc._results.lyapunov is not None:
        ly = calc._results.lyapunov
        d["lyapunov"] = {
            "V": ly.V,
            "dVdt": ly.dVdt,
            "is_stable": ly.is_stable,
            "spectral_gap": ly.spectral_gap,
        }

    return d


def kuramoto_trajectory_to_dict(trajectory: list[dict]) -> list[dict]:
    """Serialize a Kuramoto iteration trajectory for DAG inspection."""
    return [
        {
            "iteration": step["iteration"],
            "coherence": step["coherence"].tolist(),
            "K_eff": step["K_eff"].tolist(),
            "residual": step["residual"],
        }
        for step in trajectory
    ]


# ---------------------------------------------------------------------------
# Plotly figure builders
# ---------------------------------------------------------------------------

def _ensure_plotly():
    """Import plotly or raise a helpful error."""
    try:
        import plotly.graph_objects as go
        import plotly.subplots as sp
        return go, sp
    except ImportError:
        raise ImportError(
            "plotly is required for interactive visualization. "
            "Install it with: pip install plotly"
        )


def rotation_curve_figure(calc: Calculator, *,
                          show_kuramoto: bool = False,
                          title: str | None = None):
    """Interactive rotation curve with full DAG hover inspection.

    Each point shows: r, V_obs, V_MOND, a_N, a/a0, regime, coherence.
    """
    go, sp = _ensure_plotly()

    p = calc.profile
    r = p.r_kpc
    v_obs = p.v_obs
    v_mond = calc.get_rotation_curve(method="mond")
    accel = calc.get_accelerations()
    masks = calc.get_regime_mask()
    x = p.x

    regime_labels = np.array([
        "newtonian" if masks["newtonian"][i]
        else ("deep-MOND" if masks["deep_mond"][i] else "transition")
        for i in range(len(r))
    ])

    # Custom hover template showing the DAG at each radius
    hover_obs = [
        f"<b>r</b> = {r[i]:.2f} kpc<br>"
        f"<b>V_obs</b> = {v_obs[i]:.1f} km/s<br>"
        f"<b>V_MOND</b> = {v_mond[i]:.1f} km/s<br>"
        f"<b>a_N</b> = {accel['a_N'][i]:.3e} m/s²<br>"
        f"<b>a/a₀</b> = {x[i]:.3f}<br>"
        f"<b>regime</b> = {regime_labels[i]}"
        for i in range(len(r))
    ]

    fig = go.Figure()

    # Observed with error bars
    fig.add_trace(go.Scatter(
        x=r, y=v_obs,
        mode="markers",
        name="V_obs",
        marker=dict(size=6, color="#4ecdc4"),
        error_y=dict(type="data", array=p.e_v_obs, visible=True,
                     color="rgba(78,205,196,0.3)"),
        hovertext=hover_obs,
        hoverinfo="text",
    ))

    # MOND prediction — color by regime
    regime_colors = np.where(
        masks["newtonian"], "#e8e8e8",
        np.where(masks["deep_mond"], "#ff6b6b", "#ffd93d")
    )

    hover_mond = [
        f"<b>r</b> = {r[i]:.2f} kpc<br>"
        f"<b>V_MOND</b> = {v_mond[i]:.1f} km/s<br>"
        f"<b>μ(x)</b> = {accel['a_N'][i] / max(accel['a_pred'][i], 1e-30):.4f}<br>"
        f"<b>boost ν</b> = {accel['a_pred'][i] / max(accel['a_N'][i], 1e-30):.2f}×<br>"
        f"<b>deficit</b> = {max(accel['a_pred'][i] - accel['a_N'][i], 0):.3e} m/s²<br>"
        f"<b>regime</b> = {regime_labels[i]}"
        for i in range(len(r))
    ]

    fig.add_trace(go.Scatter(
        x=r, y=v_mond,
        mode="lines+markers",
        name=f"V_MOND ({calc.interpolation})",
        line=dict(color="#ff6b6b", width=2),
        marker=dict(size=4, color=regime_colors),
        hovertext=hover_mond,
        hoverinfo="text",
    ))

    # Kuramoto prediction if available
    if show_kuramoto:
        try:
            v_kur = calc.get_rotation_curve(method="kuramoto")
            ks = calc.get_kuramoto_state()
            hover_kur = [
                f"<b>r</b> = {r[i]:.2f} kpc<br>"
                f"<b>V_Kuramoto</b> = {v_kur[i]:.1f} km/s<br>"
                f"<b>coherence r(x)</b> = {ks.coherence[i]:.4f}<br>"
                f"<b>K_eff</b> = {ks.K_eff[i]:.4f}<br>"
                f"<b>ω(r)</b> = {ks.omega[i]:.4e}<br>"
                f"<b>iterations</b> = {ks.n_iter}"
                for i in range(len(r))
            ]
            fig.add_trace(go.Scatter(
                x=r, y=v_kur,
                mode="lines+markers",
                name="V_Kuramoto",
                line=dict(color="#c44dff", width=2, dash="dash"),
                marker=dict(size=4, color="#c44dff"),
                hovertext=hover_kur,
                hoverinfo="text",
            ))
        except Exception:
            pass

    fig.update_layout(
        title=title or f"Rotation Curve — {p.name}",
        xaxis_title="r [kpc]",
        yaxis_title="V [km/s]",
        template="plotly_dark",
        hovermode="closest",
        font=dict(family="monospace"),
    )

    return fig


def dag_figure(calc: Calculator, *, title: str | None = None):
    """Four-panel DAG inspection: rotation curve, RAR, Stribeck, dark matter.

    All panels share the radial axis. Hover on any panel to see the
    full computation state at that radius.
    """
    go, sp = _ensure_plotly()

    p = calc.profile
    r = p.r_kpc
    accel = calc.get_accelerations()
    x = p.x
    v_mond = calc.get_rotation_curve(method="mond")
    rho_dark = calc.get_dark_matter_density()
    prefactors = calc.get_prefactors()
    masks = calc.get_regime_mask()

    regime_labels = np.array([
        "newtonian" if masks["newtonian"][i]
        else ("deep-MOND" if masks["deep_mond"][i] else "transition")
        for i in range(len(r))
    ])

    fig = sp.make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            "Rotation Curve V(r)",
            "Radial Acceleration (RAR)",
            "Stribeck Mapping μ(x)",
            "Sync Deficit ρ_dark(r)",
        ),
        horizontal_spacing=0.12,
        vertical_spacing=0.12,
    )

    # --- Panel 1: Rotation curve ---
    fig.add_trace(go.Scatter(
        x=r, y=p.v_obs, mode="markers", name="V_obs",
        marker=dict(size=5, color="#4ecdc4"),
        error_y=dict(type="data", array=p.e_v_obs, visible=True,
                     color="rgba(78,205,196,0.3)"),
        customdata=np.column_stack([x, accel["a_N"], regime_labels]),
        hovertemplate=(
            "r=%{x:.2f} kpc<br>V_obs=%{y:.1f} km/s<br>"
            "a/a₀=%{customdata[0]:.3f}<br>"
            "a_N=%{customdata[1]:.2e}<br>"
            "regime=%{customdata[2]}<extra></extra>"
        ),
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=r, y=v_mond, mode="lines", name="V_MOND",
        line=dict(color="#ff6b6b", width=2),
    ), row=1, col=1)

    # --- Panel 2: RAR (log-log) ---
    log_aN = np.log10(np.maximum(accel["a_N"], 1e-20))
    log_aobs = np.log10(np.maximum(accel["a_obs"], 1e-20))
    log_apred = np.log10(np.maximum(accel["a_pred"], 1e-20))

    fig.add_trace(go.Scatter(
        x=log_aN, y=log_aobs, mode="markers", name="observed",
        marker=dict(size=5, color="#4ecdc4"),
        customdata=np.column_stack([r, x, prefactors]),
        hovertemplate=(
            "r=%{customdata[0]:.2f} kpc<br>"
            "log a_N=%{x:.2f}<br>log a_obs=%{y:.2f}<br>"
            "x=%{customdata[1]:.3f}<br>"
            "η=%{customdata[2]:.3f}<extra></extra>"
        ),
    ), row=1, col=2)

    fig.add_trace(go.Scatter(
        x=log_aN, y=log_apred, mode="lines", name="predicted",
        line=dict(color="#ff6b6b", width=2),
    ), row=1, col=2)

    # 1:1 line
    rar_range = [min(log_aN.min(), log_aobs.min()),
                 max(log_aN.max(), log_aobs.max())]
    fig.add_trace(go.Scatter(
        x=rar_range, y=rar_range, mode="lines", name="1:1",
        line=dict(color="gray", dash="dot", width=1),
        showlegend=False,
    ), row=1, col=2)

    # --- Panel 3: Stribeck / interpolating function ---
    x_smooth = np.logspace(-2, 2, 200)
    from sparc_x.mond import mu_rar, mu_kuramoto
    mu_rar_vals = mu_rar(x_smooth)
    mu_kur_vals = mu_kuramoto(x_smooth)
    mu_stribeck = mond_from_stribeck(x_smooth, delta=0.5)

    fig.add_trace(go.Scatter(
        x=np.log10(x_smooth), y=mu_rar_vals, mode="lines",
        name="μ_RAR", line=dict(color="#ff6b6b", width=2),
        hovertemplate="x=%{customdata:.3f}<br>μ=%{y:.4f}<extra>RAR</extra>",
        customdata=x_smooth,
    ), row=2, col=1)

    fig.add_trace(go.Scatter(
        x=np.log10(x_smooth), y=mu_kur_vals, mode="lines",
        name="μ_Kuramoto", line=dict(color="#c44dff", width=2, dash="dash"),
    ), row=2, col=1)

    fig.add_trace(go.Scatter(
        x=np.log10(x_smooth), y=mu_stribeck, mode="lines",
        name="μ_Stribeck", line=dict(color="#ffd93d", width=1, dash="dot"),
    ), row=2, col=1)

    # Overlay actual data points on the mu curve
    mu_data = accel["a_N"] / np.maximum(accel["a_pred"], 1e-30)
    fig.add_trace(go.Scatter(
        x=np.log10(np.maximum(x, 1e-10)), y=mu_data,
        mode="markers", name="data μ(x)",
        marker=dict(size=4, color="#4ecdc4"),
        customdata=np.column_stack([r, x]),
        hovertemplate=(
            "r=%{customdata[0]:.2f} kpc<br>"
            "x=%{customdata[1]:.3f}<br>"
            "μ=%{y:.4f}<extra></extra>"
        ),
    ), row=2, col=1)

    # --- Panel 4: Dark matter density (synchronization deficit) ---
    fig.add_trace(go.Scatter(
        x=r, y=np.log10(np.maximum(rho_dark, 1e-30)),
        mode="lines+markers", name="ρ_dark",
        line=dict(color="#ff6b6b", width=2),
        marker=dict(size=4),
        customdata=np.column_stack([
            rho_dark, x,
            np.maximum(accel["a_pred"] - accel["a_N"], 0),
        ]),
        hovertemplate=(
            "r=%{x:.2f} kpc<br>"
            "ρ_dark=%{customdata[0]:.3e} kg/m³<br>"
            "a/a₀=%{customdata[1]:.3f}<br>"
            "sync deficit=%{customdata[2]:.3e} m/s²"
            "<extra></extra>"
        ),
    ), row=2, col=2)

    fig.update_layout(
        title=title or f"Kuramoto DAG — {p.name}",
        template="plotly_dark",
        height=700,
        hovermode="closest",
        font=dict(family="monospace", size=11),
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=-0.15),
    )

    fig.update_xaxes(title_text="r [kpc]", row=1, col=1)
    fig.update_yaxes(title_text="V [km/s]", row=1, col=1)
    fig.update_xaxes(title_text="log₁₀ a_N", row=1, col=2)
    fig.update_yaxes(title_text="log₁₀ a_obs", row=1, col=2)
    fig.update_xaxes(title_text="log₁₀ x", row=2, col=1)
    fig.update_yaxes(title_text="μ(x)", row=2, col=1)
    fig.update_xaxes(title_text="r [kpc]", row=2, col=2)
    fig.update_yaxes(title_text="log₁₀ ρ_dark", row=2, col=2)

    return fig


def coherence_figure(calc: Calculator, *, title: str | None = None):
    """Kuramoto coherence field — the 'sheet music' of synchronization.

    Shows coherence r(x), effective coupling K_eff(x), natural frequency
    omega(x), and the Stribeck friction coefficient, all hoverable.
    """
    go, sp = _ensure_plotly()

    # Force Kuramoto solve
    ks = calc.get_kuramoto_state()
    p = calc.profile

    fig = sp.make_subplots(
        rows=2, cols=1,
        subplot_titles=(
            "Coherence r(x) — synchronization strength",
            "Effective Coupling & Natural Frequency",
        ),
        vertical_spacing=0.15,
    )

    r = p.r_kpc

    # Panel 1: coherence colored by magnitude
    fig.add_trace(go.Scatter(
        x=r, y=ks.coherence,
        mode="lines+markers",
        name="coherence r(x)",
        line=dict(color="#c44dff", width=2),
        marker=dict(
            size=6,
            color=ks.coherence,
            colorscale="Magma",
            showscale=True,
            colorbar=dict(title="r(x)", x=1.02, len=0.45, y=0.78),
        ),
        customdata=np.column_stack([
            ks.K_eff, ks.omega, p.x,
            2.0 * ks.omega,  # K_c_local
            ks.K_eff / np.maximum(2.0 * ks.omega, 1e-30),  # K/Kc ratio
        ]),
        hovertemplate=(
            "r=%{x:.2f} kpc<br>"
            "<b>coherence</b>=%{y:.4f}<br>"
            "K_eff=%{customdata[0]:.4f}<br>"
            "ω(r)=%{customdata[1]:.4e}<br>"
            "a/a₀=%{customdata[2]:.3f}<br>"
            "K_c=%{customdata[3]:.4e}<br>"
            "K/K_c=%{customdata[4]:.3f}"
            "<extra></extra>"
        ),
    ), row=1, col=1)

    # Threshold line at coherence = 0 (sync onset)
    fig.add_hline(y=0, line_dash="dot", line_color="gray",
                  annotation_text="sync onset", row=1, col=1)

    # Panel 2: K_eff and omega
    fig.add_trace(go.Scatter(
        x=r, y=ks.K_eff, mode="lines", name="K_eff",
        line=dict(color="#ffd93d", width=2),
    ), row=2, col=1)

    fig.add_trace(go.Scatter(
        x=r, y=2.0 * ks.omega, mode="lines", name="K_c = 2ω",
        line=dict(color="gray", width=1, dash="dash"),
    ), row=2, col=1)

    fig.update_layout(
        title=title or f"Coherence Field — {p.name}",
        template="plotly_dark",
        height=600,
        hovermode="closest",
        font=dict(family="monospace", size=11),
    )

    fig.update_xaxes(title_text="r [kpc]", row=1, col=1)
    fig.update_yaxes(title_text="r(x)", row=1, col=1)
    fig.update_xaxes(title_text="r [kpc]", row=2, col=1)
    fig.update_yaxes(title_text="coupling strength", row=2, col=1)

    return fig


def stribeck_figure(*, v_range: tuple[float, float] = (0.01, 10.0),
                    n_points: int = 300):
    """Standalone Stribeck curve with MOND correspondence overlay."""
    go, _ = _ensure_plotly()

    v = np.linspace(v_range[0], v_range[1], n_points)
    mu_f = stribeck_curve(v, mu_s=1.0, mu_k=0.4, v_s=1.0, delta=2.0)

    # Corresponding MOND x = v/v_s  (acceleration in units of a0)
    x_mond = v  # v_s = 1 normalisation
    mu_mond = mond_from_stribeck(x_mond, delta=0.5)

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=v, y=mu_f, mode="lines", name="Stribeck μ_f(v)",
        line=dict(color="#ff6b6b", width=2),
        customdata=np.column_stack([x_mond, mu_mond]),
        hovertemplate=(
            "v/v_s=%{x:.3f}<br>"
            "μ_friction=%{y:.4f}<br>"
            "x=a/a₀=%{customdata[0]:.3f}<br>"
            "μ_MOND=%{customdata[1]:.4f}<extra>Stribeck</extra>"
        ),
    ))

    fig.add_trace(go.Scatter(
        x=v, y=mu_mond, mode="lines", name="MOND μ(x)",
        line=dict(color="#c44dff", width=2, dash="dash"),
    ))

    fig.add_vline(x=1.0, line_dash="dot", line_color="#ffd93d",
                  annotation_text="a₀ / v_s threshold")

    fig.update_layout(
        title="Stribeck ↔ MOND Correspondence",
        xaxis_title="v / v_s  (≡ a / a₀)",
        yaxis_title="μ",
        template="plotly_dark",
        hovermode="closest",
        font=dict(family="monospace"),
    )

    return fig


# ---------------------------------------------------------------------------
# Trajectory recording for DAG animation
# ---------------------------------------------------------------------------

def solve_with_trajectory(calc: Calculator, *,
                          max_iter: int = 500,
                          tol: float = 1e-8,
                          ) -> list[dict]:
    """Re-solve Kuramoto recording every iteration for DAG inspection.

    Returns a list of dicts, one per iteration, each containing:
        iteration, coherence, K_eff, residual
    """
    from sparc_x.kuramoto import KuramotoSolver

    solver = KuramotoSolver(
        r_grid=calc.profile.r_kpc,
        omega=calc.profile.omega,
        K_coupling=calc.K_coupling,
        damping=calc.damping,
    )

    r = np.full_like(solver.r_grid, 0.5)
    alpha = solver.damping
    trajectory: list[dict] = []

    # Record initial state
    trajectory.append({
        "iteration": 0,
        "coherence": r.copy(),
        "K_eff": solver._mean_field(r).copy(),
        "residual": float("inf"),
    })

    for n in range(1, max_iter + 1):
        r_new = solver._self_consistency_step(r)
        r_update = alpha * r_new + (1.0 - alpha) * r
        residual = float(np.max(np.abs(r_update - r)))
        r = r_update

        trajectory.append({
            "iteration": n,
            "coherence": r.copy(),
            "K_eff": solver._mean_field(r).copy(),
            "residual": residual,
        })

        if residual < tol:
            break

    return trajectory


def trajectory_figure(trajectory: list[dict], r_grid: np.ndarray,
                      *, title: str | None = None):
    """Heatmap of coherence evolution — the 'score' playing out.

    x-axis: iteration (time), y-axis: radius, color: coherence r(x).
    """
    go, sp = _ensure_plotly()

    n_iter = len(trajectory)
    n_r = len(r_grid)

    # Build coherence matrix  [n_r x n_iter]
    coherence_matrix = np.zeros((n_r, n_iter))
    residuals = np.zeros(n_iter)
    for i, step in enumerate(trajectory):
        coherence_matrix[:, i] = step["coherence"]
        residuals[i] = step["residual"]

    fig = sp.make_subplots(
        rows=2, cols=1,
        subplot_titles=(
            "Synchronization Score — coherence r(x) over iterations",
            "Convergence — residual per iteration",
        ),
        row_heights=[0.75, 0.25],
        vertical_spacing=0.12,
    )

    fig.add_trace(go.Heatmap(
        z=coherence_matrix,
        x=list(range(n_iter)),
        y=r_grid,
        colorscale="Magma",
        colorbar=dict(title="r(x)", len=0.65, y=0.72),
        hovertemplate=(
            "iteration=%{x}<br>"
            "r=%{y:.2f} kpc<br>"
            "coherence=%{z:.4f}<extra></extra>"
        ),
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=list(range(n_iter)),
        y=np.log10(np.maximum(residuals, 1e-20)),
        mode="lines+markers",
        name="log₁₀ residual",
        line=dict(color="#4ecdc4", width=2),
        marker=dict(size=3),
    ), row=2, col=1)

    fig.add_hline(y=np.log10(1e-8), line_dash="dot", line_color="#ffd93d",
                  annotation_text="convergence tol", row=2, col=1)

    fig.update_layout(
        title=title or "Kuramoto DAG — Iteration Trajectory",
        template="plotly_dark",
        height=650,
        font=dict(family="monospace", size=11),
    )

    fig.update_xaxes(title_text="iteration", row=1, col=1)
    fig.update_yaxes(title_text="r [kpc]", row=1, col=1)
    fig.update_xaxes(title_text="iteration", row=2, col=1)
    fig.update_yaxes(title_text="log₁₀ |residual|", row=2, col=1)

    return fig


# ---------------------------------------------------------------------------
# Convenience: full inspection from a single call
# ---------------------------------------------------------------------------

def inspect_galaxy(calc: Calculator, *,
                   show_trajectory: bool = True,
                   export_html: str | None = None):
    """Generate all inspection figures for a galaxy.

    Parameters
    ----------
    calc : Calculator
        Configured calculator (profile already attached).
    show_trajectory : bool
        If True, solve Kuramoto and show iteration heatmap.
    export_html : str, optional
        If given, write a standalone HTML file with all figures.

    Returns
    -------
    dict of figure name -> plotly Figure
    """
    go, _ = _ensure_plotly()

    figs = {
        "rotation_curve": rotation_curve_figure(calc, show_kuramoto=show_trajectory),
        "dag": dag_figure(calc),
        "stribeck": stribeck_figure(),
    }

    if show_trajectory:
        trajectory = solve_with_trajectory(calc)
        figs["coherence"] = coherence_figure(calc)
        figs["trajectory"] = trajectory_figure(
            trajectory, calc.profile.r_kpc,
            title=f"Kuramoto DAG Trajectory — {calc.profile.name}",
        )

    if export_html:
        _export_combined_html(figs, export_html)

    return figs


def _export_combined_html(figs: dict, path: str) -> None:
    """Write all figures into a single standalone HTML file."""
    import plotly.io as pio

    parts = [
        "<!DOCTYPE html>",
        "<html><head>",
        '<meta charset="utf-8">',
        "<title>Kuramoto DAG Inspection</title>",
        "<style>body{background:#1a1a2e;color:#e8e8e8;font-family:monospace;"
        "max-width:1200px;margin:0 auto;padding:20px}"
        ".fig-container{margin:30px 0}</style>",
        "</head><body>",
        "<h1>Kuramoto DAG — Inspectable Visualizations</h1>",
    ]

    for name, fig in figs.items():
        div_html = pio.to_html(fig, full_html=False, include_plotlyjs="cdn")
        parts.append(f'<div class="fig-container"><h2>{name}</h2>{div_html}</div>')

    parts.append("</body></html>")

    with open(path, "w") as f:
        f.write("\n".join(parts))
