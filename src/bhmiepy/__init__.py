"""bhmiepy: Python interface to the BHMIE Mie scattering code.

The wrapped Fortran implementation is the BHMIE subroutine maintained by
Thomas P. Robitaille in https://github.com/hyperion-rt/bhmie (BSD-2-Clause),
itself derived from Appendix A of Bohren & Huffman (1983) with extensive
modifications by B. T. Draine. See NOTICE for attribution details.
"""

from .mie import MieResult, bhmie, compute

__version__ = "0.1.0"

__all__ = ["MieResult", "bhmie", "compute", "__version__"]
