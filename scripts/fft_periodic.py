import numpy as np

n = 10
dx = 1

def compute_fft(b):
    return np.fft.fft(b)

def compute_ifft(b):
    return np.fft.ifft(b)

def laplacian_finite_difference(b, dx):
    laplacian = np.zeros_like(b, dtype=np.float64)
    for i in range(1,n-1):
        laplacian[i] = (b[i-1] - 2*b[i] + b[i+1])/(dx**2)
    laplacian[0] = (b[-1] -2*b[0] + b[1])/(dx**2)
    laplacian[-1] = (b[-2] -2*b[-1] + b[0])/(dx**2) 
    return laplacian

# Main execution
u_ex = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])  # Example input vector
lap = laplacian_finite_difference(u_ex, dx)
print("Laplacian: ",lap)
b_hat = compute_fft(lap)
eigs = [-(2 * np.sin(k * np.pi / n) / dx)**2 for k in range(n)]
b_hat[0] = 0.0
b_hat[1:] = b_hat[1:] / eigs[1:]
b_hat = b_hat
x = compute_ifft(b_hat)

#SOLUTION IS CORRECT UP TO A CONSTANT
x = x+5.5

print("Solution: ", x.real)
