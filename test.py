from math import gcd

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

# Automatically compute factors of a given number
n = 5959  # You can change this number to any composite number
factors = pollards_rho(n)

# Print the result
print(f"Factors of {n}: {factors}")
