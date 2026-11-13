# Use Lagrange's interpolation formula to find the best estimate for y(x = 6.7) from the data table below.
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

def lagrange_interpolation(X, x):
    # Data can be either in list of (xi, yi) or just nested list (or matrices)
    if len(X) == 0:
        raise ValueError("The interpolation dataset cannot be empty.")
    
    # If in the array format
    if len(X) == 2 and not isinstance(X[0][0], (list, tuple)):
        x_values = X[0]
        y_values = X[1]
    
    # If 2 column matrix
    else:
        if any(len(row) != 2 for row in X):
            raise ValueError("Dataset must contain (x, y) pairs.")

        x_values = [row[0] for row in X]
        y_values = [row[1] for row in X]

    # Checking for length of both x and y
    if len(x_values) != len(y_values):
        raise ValueError("The x and y data must have the same length.")
    
    # Check whether the x values are distinct
    if len(set(x_values)) != len(x_values):
        raise ValueError("The x values must be distinct for Lagrange interpolation.")

    coefficients = [0.0] * len(x_values)
    for i, x_i in enumerate(x_values):
        basis = [1.0]
        for j, x_j in enumerate(x_values):
            if i != j:
                factor = [-x_j / (x_i - x_j), 1.0 / (x_i - x_j)]
                product = [0.0] * (len(basis) + 1)
                for power, coefficient in enumerate(basis):
                    product[power] += coefficient * factor[0]
                    product[power + 1] += coefficient * factor[1]
                basis = product

        for power, coefficient in enumerate(basis):
            coefficients[power] += y_values[i] * coefficient

    value = sum(coefficient * x ** power
                for power, coefficient in enumerate(coefficients))

    # Coefficients are returned in ascending order: c0, c1, ..., cn.
    return value, coefficients


# End of Code
