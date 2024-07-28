import numpy as np
import numpy.random as npr
import math
from scipy.stats import norm

def analytic_BS(sigma, r, K, T, S0):
    d_plus = 1 / (sigma * math.sqrt(T)) * (math.log(S0 / K) + (r + 0.5 * sigma * sigma) * T)
    d_minus = d_plus - sigma * math.sqrt(T)
    return S0 * norm.cdf(d_plus) - K * math.exp(-r * T) * norm.cdf(d_minus)

def european_payoff(x, K):
    return max(x - K, 0)

def generate_path_by_EM(sigma, r, K, T, S0, N):
    dt = T / N
    Z = npr.standard_normal(N)
    path = np.ones(N+1) * S0
    for n in range(1, len(path), 1):
        path[n] = path[n-1] * (1 + r * dt + sigma * math.sqrt(dt) * Z[n-1])
    return path

def montecalro(sigma, r, K, T, S0, M, N):
    exp = 0
    for m in range(M):
        exp += max(generate_path_by_EM(sigma, r, K, T, S0, N)[N] - K, 0)
    return exp / M * math.exp(-r * T)

if __name__ == "__main__":
    S0 = 125
    K = 100
    r = 0.02
    sigma = 0.3
    T = 2
    M = 10**6
    N = 100
    sim_price = monte(sigma, r, K, T, S0, M, N)
    true_price = analytic_BS(sigma, r, K, T, S0)
    error = abs(true_price - sim_price)
    print(f"analytic price : {true_price}, fdm price : {sim_price}, error : {error}")