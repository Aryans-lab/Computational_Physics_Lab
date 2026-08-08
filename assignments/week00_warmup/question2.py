"""
question2.py
------------
Problem : Calculate the sum of N terms of a Geometric Progression (GP) and
          a Harmonic Progression (HP). All parameters are read from an
          external input file.  Closed-form formulae are NOT used.
Usage   : python question2.py <input_file> <output_file>
          e.g.: python question2.py data/q2_input.txt output/q2_output.txt
Input file format (one value per line, comments with # are ignored):
    N    - number of terms
    t0   - first term
    r    - common ratio (GP)
    d    - common difference of the underlying AP (HP)
Author  : Aryan Bandyopadhyay
Roll No.: 2411014
Course  : PHY341/745 - Physics Computer Lab, NISER
"""

import sys


def sum_gp(n: int, first: float, ratio: float) -> float:
    """
    Sum the first n terms of a GP iteratively.
    t_k = first * ratio^(k-1)
    """
    total = 0.0
    term  = first
    for _ in range(n):
        total += term
        term  *= ratio
    return total


def sum_hp(n: int, first: float, diff: float) -> float:
    """
    Sum the first n terms of a HP iteratively.
    The k-th HP term is 1/a_k where a_k = first + (k-1)*diff.
    """
    total   = 0.0
    ap_term = first
    for _ in range(n):
        total   += 1.0 / ap_term
        ap_term += diff
    return total


def main() -> None:
    # -- Argument check
    if len(sys.argv) != 3:
        print("Usage: python question2.py <input_file> <output_file>",
              file=sys.stderr)
        sys.exit(1)

    input_file  = sys.argv[1]
    output_file = sys.argv[2]

    # -- Read parameters from file (skip comment lines)
    values = []
    with open(input_file, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                values.append(float(stripped.split()[0]))

    N  = int(values[0])    # number of terms
    t0 = values[1]         # first term
    r  = values[2]         # common ratio (GP)
    d  = values[3]         # common difference (HP)

    # -- Compute
    gp_sum = sum_gp(N, t0, r)
    hp_sum = sum_hp(N, t0, d)

    # -- Write output to file
    result_text = (
        f"Parameters: N={N}, t0={t0}, r={r}, d={d}\n"
        f"Sum of GP series ({N} terms): {gp_sum}\n"
        f"Sum of HP series ({N} terms): {hp_sum}\n"
    )
    with open(output_file, "w") as out:
        out.write(result_text)

    print(f"Output written to: {output_file}")


if __name__ == "__main__":
    main()

# -- Output (for N=15, t0=1.25, r=0.5, d=1.5)
# Parameters: N=15, t0=1.25, r=0.5, d=1.5
# Sum of GP series (15 terms): 2.4999237060546875
# Sum of HP series (15 terms): 2.4139570733659186
