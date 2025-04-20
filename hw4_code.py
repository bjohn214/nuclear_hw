import numpy as np
import matplotlib.pyplot as plt

def BCS(N, K=100, V=0.2, initial_mu=0.2, initial_Delta=0.2, tol=1e-6):

    def sqrt_func(epsilon_vec, Delta):
       return  np.sqrt(epsilon_vec**2 + Delta**2)

    Delta = initial_Delta
    mu = initial_mu
    k_vec = np.arange(0, K-1)
    Delta_diff = 100
    mu_diff = 100

    loop_counter = 0
    while Delta_diff >= tol and mu_diff >= tol:
        epsilon_vec = k_vec - mu
        new_Delta = (V/2)*np.sum(1/sqrt_func(epsilon_vec, Delta))
        
        new_mu = (mu/N)*np.sum(1-epsilon_vec/sqrt_func(epsilon_vec, Delta))
        Delta_diff = np.abs(Delta-new_Delta) 
        mu_diff = np.abs(mu-new_mu)

        mu = new_mu 
        Delta = new_Delta
        loop_counter += 1 


    return mu, Delta

N = np.arange(10, 42, 2)

mu = np.zeros_like(N, dtype=np.float64)
Delta = np.zeros_like(N)

for i, cN in enumerate(N): 
    mu[i] = BCS(cN)[0]
    Delta[i] = BCS(cN)[1]  

for i in range(1, len(mu)):
    if mu[i]-mu[i-1] > 0 and mu[i]-mu[i+1] > 0:
#        plt.scatter(N[i], mu[i], color='r')
        print(N[i])

plt.plot(N, mu, c='tab:blue')
plt.xlabel('N')
plt.ylabel('$\mu/\hbar\omega$')
plt.savefig('plot1.pdf')
