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
