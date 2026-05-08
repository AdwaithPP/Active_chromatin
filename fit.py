import numpy as np
import matplotlib.pyplot as plt

# =============================================
# Load data
# =============================================

exp_data = np.loadtxt("target-distribution.dat")
sim1_data = np.loadtxt("distribution-1.dat")
sim2_data = np.loadtxt("distribution-2.dat")
sim3_data = np.loadtxt("distribution-3.dat")




x = exp_data[:,0]
y = exp_data[:,1]

p1 = sim1_data[:,1]
p2 = sim2_data[:,1]
p3 = sim3_data[:,1]

# =============================================
# Select region
# =============================================

mask = x > 0

x_fit = x[mask]
y_fit = y[mask]

p1_fit = p1[mask]
p2_fit = p2[mask]
p3_fit = p3[mask]

# =============================================
# PROBABILITY-CONSTRAINED FIT (A+B+C=1)
# =============================================

best_error = 1e100
best_weights = None

a_vals = np.linspace(0, 1, 50)
b_vals = np.linspace(0, 1, 50)

for a in a_vals:
    for b in b_vals:
        c = 1 - a - b
        if c < 0:
            continue

        y_pred = a*p1_fit + b*p2_fit + c*p3_fit
        error = np.sum((y_fit - y_pred)**2)

        if error < best_error:
            best_error = error
            best_weights = (a, b, c)

A, B, C = best_weights

print("\nBest Probability Weights\n")
print("A =", A)
print("B =", B)
print("C =", C)
print("Fit error =", best_error)

