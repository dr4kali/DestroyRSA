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
        # Calculate t^2 - n
        t2_minus_n = t * t - n
        s = isqrt(t2_minus_n)

        # Check if s^2 equals t^2 - n
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

    # Handle the case for even numbers first
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

    # Loop until a non-trivial factor is found
    while d == 1:
        x = f(x)
        y = f(f(y))
        d = gcd(abs(x - y), n)

    if d == n:  # If no factor found, n is likely prime
        factors.append(n)
        return factors

    # Store the non-trivial factor found
    while n > 1:
        while n % d == 0:  # Factor out the found factor completely
            factors.append(d)
            n //= d
        if n <= 1:
            break
        # Try to find another factor
        d = pollards_rho(n)
        if d:  # Check if d has any factors returned
            d = d[0]  # Get the first factor from the returned list
        else:
            factors.append(n)  # If no factors, n is likely prime
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

        if P == Q:  # Point doubling
            if y1 == 0:  # If y is zero, return point at infinity
                return None
            m = (3 * x1**2 + self.a) * modular_inverse(2 * y1, self.p) % self.p
        else:  # Point addition
            if x1 == x2:
                return None  # Points are the same x, return point at infinity
            # Avoid division by zero
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

    # Randomly generate parameters for the elliptic curve
    for _ in range(10):  # Retry with different curves
        a = random.randint(1, n-1)
        b = random.randint(1, n-1)
        curve = EllipticCurve(a, b, n)

        # Choose a random point on the curve
        while True:
            x1 = random.randint(1, n-1)
            y1_squared = (x1**3 + a * x1 + b) % n
            y1 = isqrt(y1_squared)  # Attempt to find y
            if curve.is_point_on_curve(x1, y1):
                break

        # Start the factorization process
        P = (x1, y1)
        Q = P  # Initialize the second point
        d = 1

        for i in range(max_iterations):
            Q = curve.point_addition(Q, P)  # Q = Q + P
            if Q is None:  # If Q becomes point at infinity, break
                break
            
            d = gcd(abs(Q[0] - P[0]), n)  # GCD with the difference of x-coordinates

            if d > 1 and d not in factors:  # Found a nontrivial factor
                factors.append(d)

            # If we found all possible factors, break early
            if len(factors) >= 2:  # Adjust this threshold as needed
                break

    return factors if factors else None  # Return found factors or None

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

def general_number_field_sieve(n):
    """Simplified high-level approach of the GNFS for factoring n."""
    
    def polynomial_selection(n):
        """Step 1: Select a polynomial for the sieve."""
        # Using a simple degree-1 polynomial for educational purposes
        return lambda x: x**2 - n
    
    def sieve(polynomial, n):
        """Step 2: Perform the sieve step, collecting potential smooth relations."""
        relations = []
        B = int(math.exp(math.sqrt(math.log(n) * math.log(math.log(n)))))
        
        for x in range(1, B):
            f_x = polynomial(x)
            smooth_factors = trial_division(abs(f_x))
            if smooth_factors:
                relations.append((x, smooth_factors))
        
        return relations
    
    def trial_division(n):
        """Simple trial division to check if number is B-smooth."""
        factors = []
        for i in range(2, int(math.sqrt(n)) + 1):
            while n % i == 0:
                factors.append(i)
                n //= i
        if n > 1:
            factors.append(n)
        return factors if len(factors) > 1 else None
    
    def matrix_step(relations):
        return random.choice(relations)[0]  # Dummy relation for simplicity
    
    def square_root_step(x, n):
        """Step 4: Compute the square roots and extract factors."""
        factor = gcd(x, n)
        if factor != 1 and factor != n:
            return factor
        return None
    
    # Step 1: Polynomial Selection
    polynomial = polynomial_selection(n)
    
    # Step 2: Sieve to collect relations
    relations = sieve(polynomial, n)
    
    if not relations:
        return f"No factors found for {n}"
    
    # Step 3: Matrix step to find a dependency
    dependency = matrix_step(relations)
    
    # Step 4: Square root step to extract factors
    factors = []
    factor = square_root_step(dependency, n)
    
    if factor:
        factors.append(factor)
        factors.append(n // factor)
    
    return factors if factors else f"No factors found for {n}"

def choose_algorithm():
    """Display algorithms and let user choose."""
    algorithms = {
        '1': ('Trial Division', trial_division),
        '2': ('Fermat Factorization', fermat_factorization),
        '3': ('Pollard\'s Rho', pollards_rho),
        '4': ('Elliptic Curve Factorization', elliptic_curve_factorization),
        '5': ('Quadratic Sieve', quadratic_sieve),
        '6': ('General Number Field Sieve', general_number_field_sieve)
    }

    print("Choose a factorization algorithm:")
    for key, (name, _) in algorithms.items():
        print(f"{key}: {name}")

    choice = input("Enter the number of your choice: ")
    return algorithms.get(choice, (None, None))

def main():
    parser = argparse.ArgumentParser(description='Program to break RSA algorithm by using various factorization algorithms.')
    group = parser.add_mutually_exclusive_group(required=True)  # Make one of these options mandatory
    group.add_argument('-k', '--key', type=str, help='Path to RSA Key')
    group.add_argument('-n', '--modulus', type=int, help='Modulus (n) value')
    group.add_argument('-e', '--exponent', type=int, help='Public exponent (e) value')
    args = parser.parse_args()

    n, e = None, None  # Initialize n and e

    if args.key:
        keyPath = args.key  # Extract the filename containing Key
        try:
            with open(keyPath, "rb") as f:
                key_data = f.read()
                try:
                    # Try loading as private key
                    private_key = serialization.load_pem_private_key(
                        key_data,
                        password=None,
                        backend=default_backend()
                    )
                    n = private_key.public_key().public_numbers().n
                    e = private_key.public_key().public_numbers().e
                except ValueError:
                    # If it fails, try loading as a public key
                    public_key = serialization.load_pem_public_key(
                        key_data,
                        backend=default_backend()
                    )
                    n = public_key.public_numbers().n
                    e = public_key.public_numbers().e
        except Exception as e:
            print(f"Error while reading the key file: {e}")
            return  # Exit if key extraction fails

    if args.modulus is not None:
        n = args.modulus
    if args.exponent is not None:
        e = args.exponent
    print(f'n = {n}')
    print(f'e = {e}')
    if n is None:
        print("Failed to extract n and e from the key.")
        return  # Exit if n and e are not set

    algo_name, algo_func = choose_algorithm()

    if algo_func:
        # Determine number of cores to use (example using half the available cores)
        cores = os.cpu_count() // 2
        print(f"Using {cores} cores for factorization.")
        try:
            # Factorization using multiprocessing if needed
            with multiprocessing.Pool(processes=cores) as pool:
                factors = pool.apply(algo_func, (n,))
                if factors:
                    print("")
                    print(f"[+] Factors found: {factors}")
                else:
                    print("[+] No factors found.")
        except Exception as e:
            print(f"Error during factorization: {e}")
    else:
        print("Invalid choice. Exiting.")

if __name__ == "__main__":
    main()
