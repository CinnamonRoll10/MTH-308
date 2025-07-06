# Function to perform backward substitution (solving the upper triangular system)
def backward_sub(n, a, b):
    # Initialize the solution vector with zeros
    ans = [0] * n
    # Start solving from the last row (backward substitution)
    for i in range(n - 1, -1, -1):
        # Update the right-hand side by subtracting the values above the current row
        for j in range(i + 1, n):
            b[i] -= ans[j] * a[i][j]
        # Solve for the current variable
        ans[i] = b[i] / a[i][i]
    return ans

# Function to perform Gaussian elimination to reduce the system to row echelon form
def reduce_ref(n, a, b):
    for i in range(n - 1):  # Loop through each row (except the last one)
        # Find the row with the largest absolute value in the current column for pivoting
        mx = i
        for j in range(i + 1, n):
            if abs(a[j][i]) > abs(a[mx][i]):
                mx = j
        # Swap the rows to bring the row with the largest value to the pivot position
        a[mx], a[i] = a[i], a[mx]
        b[mx], b[i] = b[i], b[mx]

        # Perform elimination for the rows below the pivot row
        for j in range(i + 1, n):
            # Calculate the factor to eliminate the current element
            t = a[j][i] / a[i][i]
            # Subtract the pivot row's multiples from the current row to eliminate the element
            for k in range(n):
                a[j][k] -= t * a[i][k]
            b[j] -= t * b[i]

# Input the size of the system (number of variables)
n = int(input())

# Initialize the augmented matrix (a) and the vector (b)
a = []
for i in range(n):
    a.append([0] * n)  # Initialize row in matrix a
b = [0] * n  # Initialize vector b

# Read the augmented matrix (a) and vector (b) from the input
for i in range(n):
    arr = list(map(float, input().split()))  # Read the entire row
    for j in range(n):
        a[i][j] = arr[j]  # Fill the matrix a
    b[i] = arr[n]  # Fill the vector b with the last column (constant terms)

# Perform Gaussian elimination to reduce the system to row echelon form
reduce_ref(n, a, b)

# Perform backward substitution to find the solution of the system
ans = backward_sub(n, a, b)

# Print the result with one decimal place
print(" ".join(f"{num:.1f}" for num in ans))
