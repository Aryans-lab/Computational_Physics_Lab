"""
main.py
-------
Front code for PHY341/745 – Physics Computer Lab, NISER.

This script handles ALL input/output and delegates every numerical
computation to functions in mylib.py.  It does not contain any algorithm
logic itself — in compliance with the course's library/front-code rule.

Usage
-----
    python main.py <mode> <input_file> <output_file>

Modes
-----
    sum     – sum all numbers in the file
    sum2    – sum second column of a two-column (index value) file
    matadd  – add two matrices (file format: n m, then A, then B)
    matmul  – multiply two matrices (file format: n p m, then A, then B)

Author  : Aryan Bandyopadhyay
Roll No.: 2411014
"""

import sys
import mylib


def main() -> None:
    # -- Argument validation
    if len(sys.argv) != 4:
        print("Usage: python main.py <sum|sum2|matadd|matmul> "
              "<input_file> <output_file>",
              file=sys.stderr)
        sys.exit(1)

    mode        = sys.argv[1].lower()
    input_file  = sys.argv[2]
    output_file = sys.argv[3]

    # -- Dispatch to library routines
    if mode == "sum":
        tokens = mylib.read_tokens(input_file)
        values = mylib.parse_numbers(tokens)
        result = mylib.sum_numbers(values)
        text   = f"{result:.6f}\n"

    elif mode == "sum2":
        # two-column file: index  value
        values = mylib.read_second_column(input_file)
        result = mylib.sum_numbers(values)
        text   = f"{result:.6f}\n"

    elif mode == "matadd":
        tokens      = mylib.read_tokens(input_file)
        n, m, A, B  = mylib.read_matrix_add(tokens)
        C           = mylib.matrix_add(A, B, n, m)
        text        = mylib.matrix_to_string(C)

    elif mode == "matmul":
        tokens         = mylib.read_tokens(input_file)
        n, p, m, A, B  = mylib.read_matrix_multiply(tokens)
        C              = mylib.matrix_multiply(A, B, n, p, m)
        text           = mylib.matrix_to_string(C)

    else:
        print(f"Unknown mode '{mode}'. "
              "Choose one of: sum, sum2, matadd, matmul.",
              file=sys.stderr)
        sys.exit(1)

    # -- Write result to output file
    with open(output_file, "w") as f:
        f.write(text)

    print(f"Output written to: {output_file}")


if __name__ == "__main__":
    main()
