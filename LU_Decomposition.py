def is_decomposable_without_pivot(A, n):
    L = [[0] * n for _ in range(n)]
    U = [[0] * n for _ in range(n)]

    for k in range(n):
        L[k][k] = 1
        U[k][k] = A[k][k] - sum(L[k][p] * U[p][k] for p in range(k))

        if (U[k][k] == 0 and k < n - 1):
            return False

        for i in range(k + 1, n):
            sum_LU = sum(L[i][p] * U[p][k] for p in range(k))
            L[i][k] = (A[i][k] - sum_LU) / U[k][k] if U[k][k] != 0 else 0

        for j in range(k + 1, n):
            sum_LU = sum(L[k][p] * U[p][j] for p in range(k))
            U[k][j] = A[k][j] - sum_LU

    return True

def lu_decomposition(A, n):
    L = [[0] * n for _ in range(n)]
    U = [[0] * n for _ in range(n)]
    P = [[1 if i == j else 0 for j in range(n)] for i in range(n)]

    use_P = is_decomposable_without_pivot(A, n)

    for k in range(n):
        max_row = max(range(k, n), key=lambda i: abs(A[i][k]))
        if A[max_row][k] == 0:
            raise ValueError("Matrix is singular and cannot be decomposed.")

        if max_row != k and not use_P:
            A[k], A[max_row] = A[max_row], A[k]
            P[k], P[max_row] = P[max_row], P[k]
            for j in range(k):
                L[k][j], L[max_row][j] = L[max_row][j], L[k][j]

        L[k][k] = 1
        U[k][k] = A[k][k] - sum(L[k][p] * U[p][k] for p in range(k))

        for i in range(k + 1, n):
            sum_LU = sum(L[i][p] * U[p][k] for p in range(k))
            L[i][k] = (A[i][k] - sum_LU) / U[k][k] if U[k][k] != 0 else 0

        for j in range(k + 1, n):
            sum_LU = sum(L[k][p] * U[p][j] for p in range(k))
            U[k][j] = A[k][j] - sum_LU

    return use_P, L, U, P

def print_matrix(matrix, name):
    print(f"{name}:")
    for row in matrix:
        print("\t".join(f"{x:.3f}".rstrip('0').rstrip('.') if isinstance(x, float) else str(x) for x in row))

n = int(input())
A = [list(map(int, input().split())) for _ in range(n)]

done, L, U, P = lu_decomposition(A, n)

if (done):
    print("Decomposition successful without pivoting")
else:
    print("Used partial pivoting to find decomposition, permuted some rows")

print_matrix(L, "Lower Triangular Matrix L")
print_matrix(U, "Upper Triangular Matrix U")
print_matrix(P, "Permutation Matrix P")

LU = [[sum(L[i][j] * U[j][k] for j in range(n)) for k in range(n)] for i in range(n)]

# check equal
print_matrix(A, "PA Matrix (Permuted A)")
print_matrix(LU, "Product Matrix of LU")
