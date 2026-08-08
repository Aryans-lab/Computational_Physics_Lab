# read_write.py
# Minimal reference snippet: read tokens from a file, convert, write results.
# Author  : Aryan Bandyopadhyay
# Roll No.: 2411014

# -- Read
with open("input.txt", "r") as f:
    tokens = f.read().split()       # list of whitespace-separated strings

# -- Convert
numbers = [float(x) for x in tokens]

# -- Write
with open("output.txt", "w") as f:
    f.write("result here\n")
