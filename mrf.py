import math
def modified_regula_falsi(f, a, b, tol=1e-7, max_iter=100, verbose=False):
    fa = f(a)
    fb = f(b)
    if fa * fb >= 0:
        if verbose:
            print("The function must have opposite signs at a and b.")
        return None

    for i in range(max_iter):
        c = b - fb * (b - a) / (fb - fa)
        fc = f(c)
        if verbose:
            print(f"Iter {i+1}: a={a:.8f}, b={b:.8f}, c={c:.8f}, f(c)={fc:.8f}")

        if abs(fc) < tol:
            if verbose:
                print(f"Root found at x = {c:.8f}")
            return c

        # Update endpoints using modified regula falsi (Illinois)
        if fa * fc < 0:
            b, fb = c, fc
            fb /= 2  # Modification: halve the value at unchanged endpoint
        else:
            a, fa = c, fc
            fa /= 2  # Modification: halve the value at unchanged endpoint

    if verbose:
        print("Maximum iterations reached.")
    return (a + b) / 2

# Example demonstration
if __name__ == "__main__":
    def example_function(x):
        """Example: f(x) = sqrt(x) - cos(x)"""
        return x**0.5 - math.cos(x)

    root = modified_regula_falsi(example_function, 0.5, 1.0, verbose=True)
    print(f"Estimated root: {root}")
