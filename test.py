import math
from functools import reduce

def gcd(a, b):
    """Calculate the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def is_quadratic_residue(n, p):
    """Check if n is a quadratic residue mod p using Euler's criterion."""
    return pow(n, (p - 1) // 2, p) == 1

def factor_base(n, B):
    """Generate the factor base of small primes that are quadratic residues modulo n."""
    base = []
    for p in range(2, B + 1):
        if is_quadratic_residue(n, p):
            base.append(p)
    return base

def sieve(n, base, m):
    """Sieve step to find numbers whose squares are smooth."""
    smooth_numbers = []
    x_values = []
    for x in range(2, m):
        z = (x**2 - n) % n
        z_smooth = True
        z_temp = z
        factorization = []
        for p in base:
            while z_temp % p == 0:
                z_temp //= p
                factorization.append(p)
        if z_temp == 1:
            smooth_numbers.append(factorization)
            x_values.append(x)
    return smooth_numbers, x_values

def gaussian_elimination_mod2(matrix):
    """Perform Gaussian elimination on the matrix mod 2."""
    rows = len(matrix)
    cols = len(matrix[0])
    
    for col in range(cols):
        # Find a row with a 1 in the current column
        pivot_row = None
        for row in range(col, rows):
            if matrix[row][col] == 1:
                pivot_row = row
                break
        
        if pivot_row is None:
            continue
        
        # Swap rows to move pivot to the top
        matrix[col], matrix[pivot_row] = matrix[pivot_row], matrix[col]
        
        # Eliminate all other 1's in this column
        for row in range(rows):
            if row != col and matrix[row][col] == 1:
                for i in range(cols):
                    matrix[row][i] ^= matrix[col][i]
    
    # Extract dependencies (rows that sum to 0)
    dependencies = []
    for row in matrix:
        if all(v == 0 for v in row):
            dependencies.append(row)
    
    return dependencies

def quadratic_sieve(n):
    """Quadratic Sieve algorithm to find factors of n."""
    B = 50  # Factor base bound
    m = int(math.isqrt(n)) + 500  # Range for x values to gather smooth numbers

    # Step 1: Generate factor base
    base = factor_base(n, B)
    
    # Step 2: Sieve to find smooth numbers
    smooth_numbers, x_values = sieve(n, base, m)

    if len(smooth_numbers) == 0:
        print("No smooth numbers found, increase range or factor base size.")
        return None

    # Step 3: Construct matrix for Gaussian elimination mod 2
    matrix = [[0] * len(base) for _ in smooth_numbers]
    for i, factors in enumerate(smooth_numbers):
        for factor in factors:
            idx = base.index(factor)
            matrix[i][idx] += 1  # Build the matrix mod 2

    # Step 4: Solve for linear dependencies (mod 2 system)
    dependencies = gaussian_elimination_mod2(matrix)

    # Step 5: Factorization from dependencies
    factors = []
    for dep in dependencies:
        x = 1
        y = 1
        for idx, is_included in enumerate(dep):
            if is_included:
                x *= x_values[idx]
                y *= reduce(lambda a, b: a * b, [base[i] for i in smooth_numbers[idx]], 1)
        
        factor1 = gcd(x - y, n)
        factor2 = gcd(x + y, n)
        
        if factor1 != 1 and factor1 != n:
            factors.extend([factor1, factor2])

    # Return the list of factors
    return factors if factors else None

# Example usage
n = 1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139
result = quadratic_sieve(n)

if result:
    print(f"Factors of {n}: {result}")
else:
    print(f"No factors found for {n}")
