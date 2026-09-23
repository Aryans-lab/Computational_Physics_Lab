# Write library functions for Gauss-Legendre and Gauss-Laguerre quadrature for numerical integration. 
# Hardwire the table for zeros and weights for Legendre and Laguerre to any order of your 
# choice (at least n = 5). There should be only one Gauss Quadrature routine and depending 
# on what you want it should choose appropriate table for Legendre and Laguerre.
# Name: Aryan Bandyopadhyay, Roll number: 2411014

# Legendre polynomials roots and weights for n = 1 to 6

GAUSS_LEGENDRE = {
    1: ([0.0],[2.0]),
    2: ([-0.5773502691896257, 0.5773502691896257],[1.0, 1.0]),
    3: ([-0.7745966692414834, 0.0, 0.7745966692414834],
        [0.5555555555555556, 0.8888888888888889, 0.5555555555555556]),
    4: ([-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526],
        [0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538]),
    5: ([-0.9061798459386640, -0.5384693101056831, 0.0,
         0.5384693101056831, 0.9061798459386640],
        [0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
         0.4786286704993665, 0.2369268850561891]),
    6: ([-0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969,
         0.6612093864662645, 0.9324695142031521],
        [0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910,
         0.3607615730481386, 0.1713244923791704])}

# Laguerre polynomials roots and weights for n = 1 to 6

GAUSS_LAGUERRE = {
    1: ([1.0],[1.0]),   
    2: ([0.5857864376269049, 3.414213562373095],[0.8535533905932737, 0.1464466094067262]),
    3: ([0.4157745567834791, 2.294280360279041, 6.289945082937479],
        [0.7110930099291730, 0.2785177335692408, 0.0103892565015861]),
    4: ([0.3225476896193923, 1.7457611011583466, 4.536620296921128, 9.395070912301133],
        [0.6031541043416336, 0.3574186924377997, 0.0388879085150054, 0.0005392947055613]),
    5: ([0.2635603197181409, 1.413403059106517, 3.596425771040722, 7.085810005858837, 12.640800844275782],
        [0.5217556105828087, 0.3986668110831759, 0.0759424496817076, 0.0036117586799220, 0.0000233699723858]),
    6: ([0.2228466041792607, 1.188932101672623, 2.992736326059314, 5.775143569104510, 9.837467418382589, 15.982873980601701],
        [0.4589646739499636, 0.4170008307721209, 0.1133733820740449, 0.0103991974531491, 0.0002610172028149, 0.0000000000000000])}

# Gauss quadrature routine that can handle both Legendre and Laguerre methods

def gaussian_quadrature(f, a, b, N, method='legendre'):
    # Check if N is valid or not cause we are harcoding only 6 N's
    if N not in GAUSS_LEGENDRE and N not in GAUSS_LAGUERRE:
        raise ValueError("N must be one of: 1, 2, 3, 4, 5, 6.")

    if N in GAUSS_LEGENDRE:
        nodes_legendre = GAUSS_LEGENDRE[N][0]
        weights_legendre = GAUSS_LEGENDRE[N][1]
    else:
        nodes_legendre = None
        weights_legendre = None

    if N in GAUSS_LAGUERRE:
        nodes_laguerre = GAUSS_LAGUERRE[N][0]
        weights_laguerre = GAUSS_LAGUERRE[N][1]
    else:
        nodes_laguerre = None
        weights_laguerre = None

    total = 0.0

    if method == 'legendre':
        nodes = nodes_legendre
        weights = weights_legendre
        for i in range(N):
            t = nodes[i]

            # Transform node from [-1,1] to [a,b]
            x = ((b - a) / 2.0) * t + ((b + a) / 2.0)

            total = total + weights[i] * f(x)

        # Computing the integral
        integral = ((b - a) / 2.0) * total
    
    elif method == 'laguerre':
        nodes = nodes_laguerre
        weights = weights_laguerre
        for i in range(N):
            t= nodes[i]
            # No need for transformation in Laguerre as it is defined for [0, infinity)
            total = total + weights[i] * f(t)

        # For Laguerre, the integral is just the weighted sum
        integral = total

    else:
        raise ValueError("method must be 'legendre' or 'laguerre'.")

    return integral

# End of Code