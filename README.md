# DestroyRSA

A Python script to factor RSA moduli using various methods. This script supports multiple factoring algorithms and can read keys from a specified public key files.

## Features

- Factor RSA moduli using different algorithms:
  - Pollard's rho method
  - General number field sieve (GNFS)
  - Quadratic Sieve
  - Trial division
  - Fermat Factorization
  - Elliptic Curve Factorization
- Extract modulus (n) and exponent (e) directly from a public key file.
- Command-line interface for ease of use.

## Usage

1. Clone the repository:

   ```bash
   git clone https://github.com/dr4kali/destroyrsa.git
   cd destroyrsa
   ```

2. Run the python script

```
	python3 destroyrsa.py -h
```

## Examples

To factor a modulus directly:
```
python destroyrsa.py -n 5959 --method gnfs
```

To use a public key file:
```
python destroyrsa.py -k yourkey.pub --method pollard
```

