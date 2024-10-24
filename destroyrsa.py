import argparse
import os
import multiprocessing
from sympy import isprime
from sympy.ntheory import factorint
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.backends import default_backend
import math
from functools import reduce
import random

# Define the factorization algorithms
def trial_division(n):
    """Trial division factorization."""
    factors = []
    for i in range(2, int(n**0.5) + 1):
        while n % i == 0:
            factors.append(i)
            n //= i
    if n > 1:
        factors.append(n)
    return factors

def isqrt(n):
    """Integer square root using Newton's method."""
    if n < 2:
        return n
    x = n
    y = (x + n // x) // 2
    while y < x:
        x = y
        y = (x + n // x) // 2
    return x

def fermat_factorization(n):
    """Fermat's factorization method."""
    if n % 2 == 0:
        return 2, n // 2

    t = isqrt(n) + 1  # Start with the ceiling of sqrt(n)
    while True:
        t2_minus_n = t * t - n
        s = isqrt(t2_minus_n)

        if s * s == t2_minus_n:
            p = t + s
            q = t - s
            return p, q
        t += 1

def pollards_rho(n):
    """Pollard's Rho algorithm to find all factors of n."""
    if n <= 1:
        return []  # Return an empty list for non-positive integers
    
    factors = []  # List to store the factors

    while n % 2 == 0:
        factors.append(2)
        n //= 2
    
    if n == 1:
        return factors  # Return if n is fully factorized

    def f(x):
        return (x**2 + 1) % n

    x = 2
    y = 2
    d = 1

    while d == 1:
        x = f(x)
        y = f(f(y))
        d = gcd(abs(x - y), n)

    if d == n:
        factors.append(n)
        return factors

    while n > 1:
        while n % d == 0:
            factors.append(d)
            n //= d
        if n <= 1:
            break
        d = pollards_rho(n)
        if d:
            d = d[0]
        else:
            factors.append(n)
            break

    return factors

def gcd(a, b):
    """Compute GCD of a and b."""
    while b:
        a, b = b, a % b
    return a

def modular_inverse(a, m):
    """Return the modular inverse of a mod m using Extended Euclidean Algorithm."""
    m0, x0, x1 = m, 0, 1
    if m == 1:
        return 0
    while a > 1:
        q = a // m
        m, a = a % m, m
        x0, x1 = x1 - q * x0, x0
    return x1 + m0 if x1 < 0 else x1

class EllipticCurve:
    def __init__(self, a, b, p):
        self.a = a
        self.b = b
        self.p = p

    def is_point_on_curve(self, x, y):
        """Check if the point (x, y) is on the elliptic curve."""
        return (y**2 % self.p) == (x**3 + self.a * x + self.b) % self.p

    def point_addition(self, P, Q):
        """Add two points P and Q on the elliptic curve."""
        if P is None:
            return Q
        if Q is None:
            return P

        x1, y1 = P
        x2, y2 = Q

        if P == Q:
            if y1 == 0:
                return None
            m = (3 * x1**2 + self.a) * modular_inverse(2 * y1, self.p) % self.p
        else:
            if x1 == x2:
                return None
            if x2 == x1:
                return None
            m = (y2 - y1) * modular_inverse(x2 - x1, self.p) % self.p

        x3 = (m**2 - x1 - x2) % self.p
        y3 = (m * (x1 - x3) - y1) % self.p
        return (x3, y3)

def elliptic_curve_factorization(n, max_iterations=1000):
    """Elliptic Curve Factorization method to find all factors."""
    if n % 2 == 0:
        return [2]

    factors = []  # List to store found factors

    for _ in range(10):  # Retry with different curves
        a = random.randint(1, n-1)
        b = random.randint(1, n-1)
        curve = EllipticCurve(a, b, n)

        while True:
            x1 = random.randint(1, n-1)
            y1_squared = (x1**3 + a * x1 + b) % n
            y1 = isqrt(y1_squared)
            if curve.is_point_on_curve(x1, y1):
                break

        P = (x1, y1)
        Q = P
        d = 1

        for i in range(max_iterations):
            Q = curve.point_addition(Q, P)
            if Q is None:
                break
            
            d = gcd(abs(Q[0] - P[0]), n)

            if d > 1 and d not in factors:
                factors.append(d)

            if len(factors) >= 2:
                break

    return factors if factors else None

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
        pivot_row = None
        for row in range(col, rows):
            if matrix[row][col] == 1:
                pivot_row = row
                break
        
        if pivot_row is None:
            continue
        
        matrix[col], matrix[pivot_row] = matrix[pivot_row], matrix[col]
        
        for row in range(rows):
            if row != col and matrix[row][col] == 1:
                for i in range(cols):
                    matrix[row][i] ^= matrix[col][i]
    
    dependencies = []
    for row in matrix:
        if all(v == 0 for v in row):
            dependencies.append(row)
    
    return dependencies

def quadratic_sieve(n):
    """Quadratic Sieve algorithm to find factors of n."""
    B = 50  # Factor base bound
    m = int(math.isqrt(n)) + 500  # Range for x values to gather smooth numbers

    base = factor_base(n, B)
    smooth_numbers, x_values = sieve(n, base, m)

    if len(smooth_numbers) == 0:
        print("No smooth numbers found, increase range or factor base size.")
        return None

    matrix = [[0] * len(base) for _ in smooth_numbers]
    for i, factors in enumerate(smooth_numbers):
        for factor in factors:
            idx = base.index(factor)
            matrix[i][idx] += 1  # Build the matrix mod 2

    dependencies = gaussian_elimination_mod2(matrix)

    factors = []
    for dep in dependencies:
        x = 1
        y = 1
        for idx, is_included in enumerate(dep):
            if is_included:
                x *= x_values[idx]
                y *= reduce(lambda a, b: a * b, [base[i] for i in smooth_numbers[idx]], 1)
        factor = gcd(x + y, n)
        if factor != 1 and factor != n:
            factors.append(factor)

    return factors

def general_number_field_sieve(n):
    """General Number Field Sieve (GNFS) for factorization."""
    if n < 2:
        return None
    
    if isprime(n):
        return [n]  # n is prime
    
    # This is a simplified outline of the GNFS.
    # A complete implementation is complex and involves multiple steps.
    
    # Step 1: Select a polynomial.
    # (We would select two polynomials; for simplicity, we'll use one.)
    f = lambda x: x**2 - n

    # Step 2: Find a smooth number relation.
    # In practice, this involves sieving.
    smooth_numbers = []  # Placeholder for smooth number relations

    # Simulating smooth numbers for the example
    for x in range(2, 5000):
        z = f(x)
        if isprime(z):  # Simulate finding smooth numbers
            smooth_numbers.append(z)

    if len(smooth_numbers) < 2:
        print("Not enough smooth numbers found.")
        return None

    # Step 3: Use linear algebra over GF(2) to find dependencies.
    matrix = [[0] * len(smooth_numbers) for _ in range(len(smooth_numbers))]
    for i in range(len(smooth_numbers)):
        for j in range(len(smooth_numbers)):
            # Placeholders for linear dependencies.
            matrix[i][j] = random.randint(0, 1)

    # Step 4: Perform Gaussian elimination to find dependencies.
    dependencies = gaussian_elimination_mod2(matrix)

    factors = []
    for dep in dependencies:
        # Placeholder to calculate factors from dependencies.
        factor = reduce(lambda x, y: x * y, [smooth_numbers[i] for i in range(len(dep)) if dep[i]], 1)
        if factor not in factors and factor != 1 and factor != n:
            factors.append(factor)

    return factors

def extract_key_data(public_key_file):
    """Extract modulus n and exponent e from a public key file."""
    if not os.path.isfile(public_key_file):
        raise FileNotFoundError(f"The file {public_key_file} does not exist.")

    with open(public_key_file, 'rb') as key_file:
        public_key = serialization.load_pem_public_key(
            key_file.read(),
            backend=default_backend()
        )

    # Extracting the modulus n and exponent e
    numbers = public_key.public_numbers()
    n = numbers.n
    e = numbers.e
    return n, e

def factor(n, method='trial'):
    """Main function to factor n using various algorithms."""
    if n < 2:
        return []

    # Try using the selected method
    if method == 'trial':
        return trial_division(n)
    elif method == 'fermat':
        return fermat_factorization(n)
    elif method == 'pollards_rho':
        return pollards_rho(n)
    elif method == 'elliptic_curve':
        return elliptic_curve_factorization(n)
    elif method == 'quadratic_sieve':
        return quadratic_sieve(n)
    elif method == 'gnfs':
        return general_number_field_sieve(n)
    else:
        raise ValueError("Unknown method specified.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Factor an integer.')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-k', '--key', type=str, help='Public key file to extract n and e.')
    group.add_argument('-n', '--modulus', type=int, help='The modulus n to factor.')
    parser.add_argument('-e', '--exponent', type=int, help='The exponent e (optional).')
    parser.add_argument('--method', type=str, choices=['trial', 'fermat', 'pollards_rho', 'elliptic_curve', 'quadratic_sieve', 'gnfs'], default='trial', help='Factorization method to use.')

    args = parser.parse_args()

    # Determine n from key or modulus
    if args.key:
        n, e = extract_key_data(args.key)
    elif args.modulus is not None:
        n = args.modulus
    else:
        raise ValueError("You must specify either a public key file or modulus.")

    result = factor(n, method=args.method)
    print(f'Factors of {n} using {args.method}: {result}')
