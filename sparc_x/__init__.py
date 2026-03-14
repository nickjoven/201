"""
SPARC-X-API: Python API for the Kuramoto-Einstein synchronization framework.

Provides computational tools for:
- MOND interpolation via Stribeck friction correspondence
- Kuramoto order-parameter dynamics on a radial manifold
- SPARC galaxy rotation-curve fitting and prediction
- Lyapunov stability analysis for synchronization fixed points
"""

from sparc_x.constants import A0, G
from sparc_x.calculator import Calculator

__version__ = "0.1.0"
__all__ = ["Calculator", "A0", "G"]


def inspect(name_or_path: str = "NGC6503", **kwargs):
    """Quick-start: load a galaxy and return inspectable figures.

    >>> import sparc_x
    >>> figs = sparc_x.inspect("NGC6503")
    >>> figs["dag"].show()
    """
    from sparc_x.viz import inspect_galaxy
    from sparc_x.profiles import load_rotmod, load_sparc_table
    from pathlib import Path

    data_dir = Path(__file__).parent.parent / "data"
    dat = data_dir / f"{name_or_path}_rotmod.dat"
    if dat.exists():
        profile = load_rotmod(dat)
    else:
        table = load_sparc_table(data_dir)
        if name_or_path in table:
            profile = table[name_or_path]
        else:
            raise FileNotFoundError(
                f"Galaxy {name_or_path!r} not found. "
                f"Available: {sorted(table.keys())[:10]}..."
            )
    calc = Calculator(profile)
    return inspect_galaxy(calc, **kwargs)
