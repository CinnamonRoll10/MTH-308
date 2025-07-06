import math

# Newton-Raphson Method for solving nonlinear equations
# This implementation solves a system of nonlinear equations iteratively
# using the Jacobian matrix and Gaussian elimination.

# Function to perform backward substitution
# Solves the upper triangular system after Gaussian elimination
def backward_sub(n, a, b):
    ans = [0] * n

    for i in range(n - 1, -1, -1):
        for j in range(i + 1, n):
            b[i] -= ans[j] * a[i][j]  # Updating b[i] based on known values
        ans[i] = b[i] / a[i][i]  # Compute the solution for x_i
    return ans

# Function to convert a matrix into row echelon form
def reduce_ref(n, a, b):
    for i in range(n - 1):
        mx = i  # Index of the maximum element in column i
        for j in range(i + 1, n):
            if abs(a[j][i]) > abs(a[mx][i]):
                mx = j  # Swap rows if needed

        a[mx], a[i] = a[i], a[mx]  # Swap the maximum row with the current row
        b[mx], b[i] = b[i], b[mx]

        for j in range(i + 1, n):
            t = a[j][i] / a[i][i]  # Compute the factor for row reduction
            for k in range(n):
                a[j][k] -= t * a[i][k]  # Reduce row j
            b[j] -= t * b[i]  # Adjust the right-hand side vector

# Function to perform Gaussian elimination
def gauss_elim(n, a, b):
    reduce_ref(n, a, b)  # Convert to row echelon form
    ans = backward_sub(n, a, b)  # Solve the system
    return ans

# Function defining the system of equations
def fi(i, x):
    y = []
    y.append(0)  # Dummy value to make indexing easier
    for k in range(n):
        y.append(x[k])

    if i == 0:
        return y[1] * y[1] + y[1] * y[2] - 10  # First equation
    elif i == 1:
        return y[2] + 3 * y[1] * y[2] * y[2] - 57  # Second equation

# Function to compute partial derivatives of fi w.r.t xj
def dfidxj(i, j, x):
    y = []
    y.append(0)  # Dummy value for easier indexing
    for k in range(n):
        y.append(x[k])

    if i == 0:
        if j == 0:
            return 2 * y[1] + y[2]  # Partial derivative of first equation
        elif j == 1:
            return y[1]  # Partial derivative of first equation
    elif i == 1:
        if j == 0:
            return 3 * y[2] * y[2]  # Partial derivative of second equation
        elif j == 1:
            return 1 + 6 * y[1] * y[2]  # Partial derivative of second equation

# Taking user input for the system size
print("Enter n")
n = int(input())

# Taking initial guesses from user
print("Enter the initial guesses")
x = list(map(float, input().split()))

# Taking the required error threshold
print("Enter the needed error")
es = float(input())

it = 1  # Iteration counter

# Print the table header
print(f"{'Iteration':<15}{'x_in':<15}{'K':<35}{'f':<15}{'delta_x':<15}{'X_out':<15}{'error':<10}")

# Iterative process using the Newton-Raphson method
while(it < 30):
    # Initialize Jacobian matrices
    K = []
    K2 = []
    for i in range(n):
        empty_arr = [0] * n
        K.append(empty_arr)
        empty_arr2 = [0] * n
        K2.append(empty_arr2)

    # Compute Jacobian matrix
    for i in range(n):
        for j in range(n):
            K[i][j] = dfidxj(i, j, x)
            K2[i][j] = dfidxj(i, j, x)

    # Compute function values
    f = [0] * n
    F = [0] * n
    for i in range(n):
        F[i] = -fi(i, x)
        f[i] = -F[i]

    # Solve the linear system for delta_x
    delta = gauss_elim(n, K, F)

    # Compute new approximations
    new_x = [0] * n
    for i in range(n):
        new_x[i] = x[i] + delta[i]

    # Compute error
    wrong_sum = 0
    x_sum = 0

    for i in range(n):
        wrong_sum += delta[i] * delta[i]
        x_sum += new_x[i] * new_x[i]

    error = wrong_sum / x_sum
    error = math.sqrt(error)

    # Print the iteration details
    for i in range(n):
        x_str = format(x[i], '.5f')
        f_str = format(f[i], '.5f')
        delta_str = format(delta[i], '.5f')
        new_x_str = format(new_x[i], '.5f')
        K2_str = str([format(num, '.5f') for num in K2[i]])

        if i == 0:
            print(f"{it:<15}{str([x_str]):<15}{K2_str:<35}{str([f_str]):<15}{str([delta_str]):<15}{str([new_x_str]):<15}{format(error, '.5f'):<10}")
        else:
            print(f"{'':<15}{str([x_str]):<15}{K2_str:<35}{str([f_str]):<15}{str([delta_str]):<15}{str([new_x_str]):<15}")

    # Update the solution for the next iteration
    x = new_x
    it = it + 1

    # Check for convergence
    if (error < es):
        break

# Print the final results
print("Final results")
print([format(num, ".5f") for num in x])
