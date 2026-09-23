# Find the root of the following function using Fixed-point method f(x)=x^2-2x-3
# Name: Aryan Bandyopadhyay, Roll number: 2411014

from mylib import fixed_point

if __name__ == '__main__':
    # Equation: f(x) = x^2 - 2x - 3 = 0
    # Analytical roots are x = -1 and x = 3

    print("Fixed Point Method:")
    # For root x = 3, choose g(x) = sqrt(2x + 3)
    def g1(x):
        return (2*x + 3)**0.5
    print ("\nFixed Point Method for first root (x = 3):")
    initial_guess_fp1 = 1
    root_fp1 = fixed_point(g1, x0=initial_guess_fp1, accuracy=1e-6)

# End of Code

############################################
# Output:
# Fixed Point Method:
#
# Fixed Point Method for first root (x = 3):
# The root is: 2.9999998290170797
# The number of iterations taken is 15
##############################################