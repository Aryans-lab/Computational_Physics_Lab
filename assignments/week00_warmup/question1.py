"""
question1.py
------------
Problem : Compute the sum of the first N even numbers and the factorial of M,
          where N and M are read from an external input file.
          Uses an explicit loop -- no math.factorial or built-in sum().
Usage   : python question1.py <input_file> <output_file>
          e.g.: python question1.py data/q1_input.txt output/q1_output.txt
Author  : Aryan Bandyopadhyay
Roll No.: 2411014
Course  : PHY341/745 - Physics Computer Lab, NISER
"""

import sys


def sum_of_first_n_even_numbers(n: int) -> int:
    """Return the sum 2 + 4 + 6 + ... + 2n using a loop."""
    total = 0
    for i in range(1, n + 1):
        total += 2 * i
    return total


def factorial(n: int) -> int:
    """Return n! using a loop (no math.factorial)."""
    result = 1
    for i in range(1, n + 1):
        result *= i
    return result


def main() -> None:
    # -- Argument check
    if len(sys.argv) != 3:
        print("Usage: python question1.py <input_file> <output_file>",
              file=sys.stderr)
        sys.exit(1)

    input_file  = sys.argv[1]
    output_file = sys.argv[2]

    # -- Read N and M from file (skip comment lines starting with #)
    values = []
    with open(input_file, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                for token in stripped.split():
                    values.append(int(token))

    n_even = values[0]   # sum of first n_even even numbers
    m_fact = values[1]   # compute m_fact!

    # -- Compute
    s = sum_of_first_n_even_numbers(n_even)
    f = factorial(m_fact)

    # -- Write output to file
    result_text = (
        f"Sum of first {n_even} even numbers: {s}\n"
        f"Factorial of {m_fact}             : {f}\n"
    )
    with open(output_file, "w") as out:
        out.write(result_text)

    print(f"Output written to: {output_file}")


if __name__ == "__main__":
    main()

# -- Output (for input N=20, M=8)
# Sum of first 20 even numbers: 420
# Factorial of 8             : 40320
