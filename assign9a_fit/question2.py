# Write a library function for linear regression, which will also calculate Pearson's correlation coefficient r^2
# Name: Aryan Bandyopadhyay, Roll Number: 2411014

def linear_regression(X, sigma=None):

    if not X:
        raise ValueError("The dataset cannot be empty.")

    if len(X) == 2 and not isinstance(X[0][0], (list, tuple)):
        x_values = X[0]
        y_values = X[1]
    else:
        if any(len(row) != 2 for row in X):
            raise ValueError("Dataset must contain (x, y) pairs.")
        x_values = [row[0] for row in X]
        y_values = [row[1] for row in X]

    if len(x_values) != len(y_values):
        raise ValueError("The x and y data must have the same length.")
    if len(x_values) < 2:
        raise ValueError("At least two data points are required.")

    n = len(x_values)
    if sigma is None:
        sigma_values = [1.0] * n
    elif isinstance(sigma, (int, float)):
        sigma_values = [float(sigma)] * n
    else:
        sigma_values = list(sigma)
        if len(sigma_values) != n:
            raise ValueError("sigma must have one value for each data point.")

    if any(value <= 0 for value in sigma_values):
        raise ValueError("All sigma values must be positive.")

    weights = [1.0 / value ** 2 for value in sigma_values]
    S   = sum(weights)
    Sx  = sum(weight * xi         for weight, xi         in zip(weights, x_values))
    Sy  = sum(weight * yi         for weight, yi         in zip(weights, y_values))
    Sxx = sum(weight * xi ** 2    for weight, xi         in zip(weights, x_values))
    Sxy = sum(weight * xi * yi    for weight, xi, yi     in zip(weights, x_values, y_values))
    Syy = sum(weight * yi ** 2    for weight, yi         in zip(weights, y_values))  # FIX: was x_values

    # Calculate slope (m) and intercept (b)
    denominator = S * Sxx - Sx ** 2
    if denominator == 0:
        raise ValueError("Denominator in slope calculation is zero. Check your data points.")

    slope     = (S * Sxy - Sx * Sy) / denominator
    intercept = (Sxx * Sy - Sx * Sxy) / denominator

    # Pearson correlation coefficient r^2
    # r = (S*Sxy - Sx*Sy) / sqrt((S*Sxx - Sx^2) * (S*Syy - Sy^2))
    # r^2 = (S*Sxy - Sx*Sy)^2 / ((S*Sxx - Sx^2) * (S*Syy - Sy^2))
    r_squared = (S * Sxy - Sx * Sy) ** 2 / (denominator * (S * Syy - Sy ** 2))

    return slope, intercept, r_squared

# End of Code
