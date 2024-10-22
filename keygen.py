from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives.asymmetric import rsa
from cryptography.hazmat.primitives import serialization

# Generate a 64-bit RSA key pair
private_key = rsa.generate_private_key(
    public_exponent=65537,
    key_size=64,  # Small key size for testing
    backend=default_backend()
)

# Serialize the private key to PEM format
pem_private = private_key.private_bytes(
    encoding=serialization.Encoding.PEM,
    format=serialization.PrivateFormat.TraditionalOpenSSL,
    encryption_algorithm=serialization.NoEncryption()
)

# Serialize the public key to PEM format
public_key = private_key.public_key()
pem_public = public_key.public_bytes(
    encoding=serialization.Encoding.PEM,
    format=serialization.PublicFormat.SubjectPublicKeyInfo
)

# Save the keys to files
with open('private_key_64.pem', 'wb') as f:
    f.write(pem_private)

with open('public_key_64.pem', 'wb') as f:
    f.write(pem_public)

print("64-bit RSA keys generated and saved as private_key_64.pem and public_key_64.pem.")
