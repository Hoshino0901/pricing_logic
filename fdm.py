import numpy as np
import math
from scipy.stats import norm

def analytic_BS(sigma, r, K, T, S0):
    d_plus = 1 / (sigma * math.sqrt(T)) * (math.log(S0 / K) + (r + 0.5 * sigma * sigma) * T)
    d_minus = d_plus - sigma * math.sqrt(T)
    return S0 * norm.cdf(d_plus) - K * math.exp(-r * T) * norm.cdf(d_minus)

def terminal_condition(K, x):
    return max(x - K, 0)

def boundary_condition_at_zero(K, t):
    return 0

def boundary_condition_at_infty(r, T, X, K, t):
    return X - math.exp(-r * (T - t)) * K

def implicit_fdm(sigma, r, K, T, S0, N, x=8):
    X = S0 * x
    while True:
        if X > K:
            break
        else:
            x += 1
            X = S0 * x
    M = x * 100
    dx = X / M
    dt = T / N

    f_sim = np.zeros((N+1, M+1))
    for j in range(M+1):
        f_sim[N, j] = terminal_condition(K, j * dx)
    for i in range(N+1):
        f_sim[i, 0] = boundary_condition_at_zero(K, i * dt)
        f_sim[i, M] = boundary_condition_at_infty(r, T, X, K, i * dt)

    a = np.zeros(M) #a[0]は使わない
    b = np.zeros(M) #b[0]は使わない
    c = np.zeros(M) #c[0]は使わない
    d = np.zeros(M) #d[0]は使わない

    for i in range(N-1, -1, -1):
        for j in range(1, M):
            a[j] = -0.5 * dt * (sigma**2 * j**2 + r * j)
            b[j] = 1 + dt * (sigma**2 * j**2 + r)
            c[j] = -0.5 * dt * (sigma**2 * j**2 - r * j)
            d[j] = f_sim[i + 1, j]

        d[1] -= c[1] * f_sim[i, 0]
        d[M-1] -= a[M-1] * f_sim[i, M]

        P = np.zeros(M) #P[0]は使わない
        Q = np.zeros(M) #Q[0]は使わない
        P[1] = a[1] / b[1]
        Q[1] = d[1] / b[1]
        for j in range(2, M):
            P[j] = a[j] / (b[j] - c[j] * P[j-1])
            Q[j] = (d[j] - c[j] * Q[j-1]) / (b[j] - c[j] * P[j-1])

        f_sim[i, M-1] = Q[M-1]
        for j in range(M-2, 0, -1):
            f_sim[i, j] = -P[j] * f_sim[i, j+1] + Q[j] 

    index = int(S0 / dx)
    return f_sim[0, index]

if __name__ == "__main__":
    S0 = 125
    K = 100
    r = 0.02
    sigma = 0.3
    T = 2
    N = 5000
    sim_price = implicit_fdm(sigma, r, K, T, S0)
    true_price = analytic_BS(sigma, r, K, T, S0)
    error = abs(true_price - sim_price)
    print(f"analytic price : {true_price}, fdm price : {sim_price}, error : {error}")