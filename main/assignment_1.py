# JP O'Toole
# COT 4500
# Programming Assignment 2

import numpy as np

# Question 1
def neville1(x_points, y_points, x):

    # Records number of points
    n = len(x_points)
    p = [[0] * n for _ in range(n)]

    # Appends X and Y
    for i in range(n):
        p[i][0] = y_points[i]
    
    # Loops and does computations
    for j in range(1, n):
        for i in range(j, n):
            p[i][j] = ((x - x_points[i-j]) * p[i][j-1] - (x - x_points[i]) * p[i-1][j-1]) / (x_points[i] - x_points[i-j])

    # Returns final value
    return p[n-1][n-1]

# Values in question
x_points = [3.6, 3.8, 3.9]
y_points = [1.675, 1.436, 1.318]
x = 3.7

# Calls functions and prints value
interp = neville1(x_points, y_points, x)
print(interp, "\n")

# Question 2 & 3
def newtonfwdcoeffs(x_points, y_points):

    n = len(x_points)
    fwd_diff = np.zeros((n, n))

    # Fill first two columns with x and f(x)
    fwd_diff[:, 0] = y_points

    # Compute forward differences
    for j in range(1, n):
        for i in range(1, j+1):
            fwd_diff[j][i] = (fwd_diff[j][i-1] - fwd_diff[j-1][i-1])
            if(i == 1):
                h = x_points[j] - x_points[j-1]
                fwd_diff[j][i] = fwd_diff[j][i] / h
            elif(i == 2):
                h = x_points[j] - x_points[j-2]
                fwd_diff[j][i] = fwd_diff[j][i] / h
            else:
                h = x_points[j] - x_points[j-3]
                fwd_diff[j][i] = fwd_diff[j][i] / h
               
    print(fwd_diff)

    # Returns the polynomial coefficients
    coeffs = [fwd_diff[1,1], fwd_diff[2,2], fwd_diff[3,3]]
    return coeffs

def newtonfwdinterp(x_points, y_points, x_approx):

    # Re-establishes variables
    n = len(x_points)
    h = x_points[1] - x_points[0]
    u = (x_approx - x_points[0]) / h
    term = 1
    y_approx = y_points[0]

    # Calls coefficients function
    poly = newtonfwdcoeffs(x_points, y_points)

    # Approximates y
    for i in range(len(poly)):
        term *= (u-1) / (i+1)
        y_approx += poly[i] * term

    # Returns key values
    return poly, y_approx


# Values in Question
x_points = [7.2, 7.4, 7.5, 7.6]
y_points = [23.5492, 25.3913, 26.8224, 27.4589]
x_approx = 7.3

# Calls function
poly, y_approx = newtonfwdinterp(x_points, y_points, x_approx)

# Prints the Coefficients and Approximates for 7.3
print(poly[0], "\n", poly[1], "\n", poly[2], "\n")

# Prints f(7.3)
print(y_approx, "\n")
    
# Question 4
def hermite_divided_difference(x_vals, y_vals, dy_vals):
    n = len(x_vals)
    
    # Create Hermite matrix with duplicated x values
    Q = np.zeros((n, n+1))
    Q[:, 0] = x_vals
    Q[:, 1] = y_vals

    # First order derivatives
    for i in range(1, n):
        if x_vals[i] == x_vals[i-1]:
            Q[i, 2] = dy_vals[i]
        else:
            Q[i, 2] = (Q[i, 1] - Q[i - 1, 1]) / (x_vals[i] - x_vals[i-1])
    
    # Higher-order divided differences
    for j in range(3, n+1):
        for i in range(j-1, n):
            Q[i, j] = (Q[i, j-1] - Q[i-1, j-1]) / (x_vals[i] - x_vals[i-j+1])

    # Display the divided difference table
    for row in Q:
        print(["{:.8e}".format(val) for val in row])
    print("\n")

# Given data
x_values = [3.6, 3.6, 3.8, 3.8, 3.9, 3.9]
y_values = [1.675, 1.675, 1.436, 1.436, 1.318, 1.318]
dy_values = [-1.195, -1.195, -1.188, -1.188, -1.182, -1.182]

# Compute Hermite divided difference table
hermite_divided_difference(x_values, y_values, dy_values)

# Question 5
def cspline(x, y):

    # Establishes variables
    n = len(x)
    h = np.diff(x)
    A = np.zeros((n, n))
    b = np.zeros(n)

    # Establishes [1, 0, 0, 0] for row 1 and [0, 0, 0, 1] for row 4
    A[0, 0] = 1
    A[-1, -1] = 1

    # Performs the interpolation for Matrix A
    for i in range(1, n-1):
        A[i, i-1] = h[i-1]
        A[i, i] = 2 * (h[i-1] + h[i])
        A[i, i+1] = h[i]

    # Prints matrix A
    for row in A:
        print(row)

    # Solves for vector b and prints it
    for i in range(1, n-1):
        b[i] = ((y[i+1] - y[i]) / h[i]) - ((y[i] - y[i-1]) / h[i-1])
    print(b)

    # Solves for vector x and prints it
    x = np.linalg.solve(A, b)
    print(x)

# Sets array and calls the function
x = np.array([2, 5, 8, 10])
y = np.array([3, 5, 7, 9])
cspline(x, y)
