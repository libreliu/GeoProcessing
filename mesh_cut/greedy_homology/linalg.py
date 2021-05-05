import numpy as np

def check_z2(A: np.ndarray):
    """Check if given matrix is in Z_2"""
    return ((A == 1) + (A == 0)).all()

# https://math.stackexchange.com/questions/3073083/how-to-reduce-matrix-into-row-echelon-form-in-numpy/3073117
def get_Bopt_column(A_input: np.ndarray):
    """Get first rank(A) linear independent column vectors of A over Z_2"""
    A = A_input.copy()
    assert(check_z2(A))

    m, n = A.shape
    pivot_column = []

    # no need to care things upside working_row
    working_row = 0
    working_col = 0

    while working_row < m and working_col < n:
        # find element with first non-zero coeff this col
        for i in range(working_row, m):
            if A[i, working_col] != 0:
                break
        else:
            # all elements in this column is zero
            working_col += 1
            continue

        if i > working_row:
            # swap row working_row with row i (working row)
            A[[working_row, i], :] = A[[i, working_row], :]
        
        # (rows, cols) = (rows, cols) - row_vec * coeff_in_working_col
        for r in range(working_row + 1, m):
            A[r] ^= A[working_row] & A[r, working_col]

        pivot_column.append(working_col)
        working_row += 1
        working_col += 1

    return pivot_column

def solve_z2_sequential(A_input: np.ndarray, Z: np.ndarray):
    m, n = A_input.shape
    assert(m > n)
    row_pivot = get_Bopt_column(A_input.T)
    assert(len(row_pivot) == n)

    X = np.linalg.inv(A_input[row_pivot, :]) @ Z[row_pivot, :]

    result = (np.round(X) % 2).astype(np.int8)
    #assert(np.allclose(result, solve_z2_sequential_slow(A_input, Z)))

    return result

# TODO: pre-calculate factorization
def solve_z2_sequential_slow(A_input: np.ndarray, Z: np.ndarray):
    m, n = A_input.shape
    m_dup, p = Z.shape
    assert(m_dup == m)
    solution = np.ndarray((n, p), dtype=np.int8)
    for i in range(0, p):
        solution[:, [i]] = solve_z2(A_input, Z[:, [i]])
    
    return solution

def solve_z2(A_input: np.ndarray, z: np.ndarray):
    """solve over Z_2 coeff, z must be of dim 2
    When num of solutions > 1, only Bopt pivot column have coeff"""
    m, n = A_input.shape
    A_aug = np.hstack((A_input, z))

    #assert(n not in get_Bopt_column(A_aug))

    pivot_column = []

    # no need to care things upside working_row
    working_row = 0
    working_col = 0

    while working_row < m and working_col < n:
        # find element with first non-zero coeff this col
        for i in range(working_row, m):
            if A_aug[i, working_col] != 0:
                break
        else:
            # all elements in this column is zero
            working_col += 1
            continue

        if i > working_row:
            # swap row working_row with row i (working row)
            A_aug[[working_row, i], :] = A_aug[[i, working_row], :]

        # (rows, cols) = (rows, cols) - row_vec * coeff_in_working_col
        for r in range(working_row + 1, m):
            A_aug[r] ^= A_aug[working_row] & A_aug[r, working_col]

        for r in range(0, working_row):
            A_aug[r] ^= A_aug[working_row] & A_aug[r, working_col]

        pivot_column.append(working_col)
        working_row += 1
        working_col += 1

    # Backsubstitution from upper triangular
    assert(len(pivot_column) == n)
    solution = np.ndarray(shape=(n,1), dtype=np.int8)
    #assert((A_aug[0:n, 0:n] == np.identity(n, dtype=np.int8)).all())
    for i in range(n-1, -1, -1):
        solution[i] = A_aug[i, n]

    return solution


# def solve_lstsq(A_input: np.ndarray, b: np.ndarray):
#     m, n = A_input.shape

#     A_aug = GenericMatrix(
#             size=(n, n),
#             zeroElement=0,
#             identityElement=1,
#             add=lambda x, y: x ^ y,  # XOR
#             mul=lambda x, y: x & y,  # AND
#             sub=lambda x, y: x ^ y,  # XOR
#             div=lambda x, y: x       # Divide by 1 yields identity
#         )
    

#     for i in range(0, n):
#         for j in range(0, n):
#             for k in range(0, m):
#                 A_aug[i, j] ^= (int(A_input[k, i]) & int(A_input[k, j]))
    

#     b_aug = GenericMatrix(
#             size=(n, 1),
#             zeroElement=0,
#             identityElement=1,
#             add=lambda x, y: x ^ y,  # XOR
#             mul=lambda x, y: x & y,  # AND
#             sub=lambda x, y: x ^ y,  # XOR
#             div=lambda x, y: x       # Divide by 1 yields identity
#         )

#     for i in range(0, n):
#         for k in range(0, m):
#             b_aug[i, 0] ^= (int(A_input[k, i]) & int(b[k]))


#     x = A_aug.Solve(b_aug.GetColumn(0))
    
#     return x