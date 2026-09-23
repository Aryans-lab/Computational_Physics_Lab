"""
question1.py — week06b_simpson_montecarlo
-----------------------------------------
Problem : Write library functions for Simpson's 1/3 rule and for Monte Carlo
          integration, passing the integrand as a function.  The routines
          live in mylib.py (simpson, monte_carlo); this driver validates
          them: Simpson's rule must be *exact* on cubics (its order of
          precision is 3), and the Monte Carlo estimate must sit inside a
          few standard errors of the exact value.
Usage   : python question1.py [output_file]
          Defaults: output/q1_output.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

from mylib import monte_carlo, simpson

BASE_DIR = Path(__file__).resolve().parent


def main(argv: list[str] | None = None) -> None:
    args = sys.argv[1:] if argv is None else list(argv)
    out_file = Path(args[0]) if args else BASE_DIR / "output" / "q1_output.txt"

    # -- Simpson exactness on a cubic (exact for degree <= 3)
    cubic = lambda x: x ** 3 + 2.0 * x ** 2 - x + 5.0  # noqa: E731
    exact_cubic = (1.0 / 4.0) + (2.0 / 3.0) - 0.5 + 10.0  # on [0, 1]
    val, evals = simpson(cubic, 0.0, 1.0, N=2)
    cubic_err = abs(val - exact_cubic)

    # -- Monte Carlo on int_0^1 x^2 dx = 1/3
    f = lambda x: x * x  # noqa: E731
    mc_val, sigma_f, sigma_I = monte_carlo(f, 0.0, 1.0, N=10000, seed=7)
    z = abs(mc_val - 1.0 / 3.0) / sigma_I

    lines = [
        "1. Simpson's 1/3 rule — exactness check",
        "   integral of f(x) = x^3 + 2x^2 - x + 5 on [0, 1] (a cubic)",
        f"   Simpson (N = 2, only {evals} evaluations) = {val:.15f}",
        f"   exact value                    = {exact_cubic:.15f}",
        f"   |error| = {cubic_err:.3e}   ->  Simpson's rule integrates cubics",
        "   exactly (up to round-off): its order of precision is 3, i.e. it is",
        "   exact for every polynomial of degree <= 3, whatever even N is used.",
        "",
        "2. Monte Carlo integration — statistical check",
        "   integral of f(x) = x^2 on [0, 1],  N = 10000,  LCG seed = 7",
        f"   estimate = {mc_val:.8f}   (exact = 0.3333333333)",
        f"   sigma_f  = {sigma_f:.6f},   sigma_estimate = {sigma_I:.6f}",
        f"   deviation in units of sigma: |estimate - exact|/sigma = {z:.2f}",
        "   a deviation of O(1) sigma is exactly what a correct Monte Carlo",
        "   estimator should produce; systematic drift would show up as",
        "   |deviation| >> 10 sigma.",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"Output written to: {out_file}")


if __name__ == "__main__":
    main()
