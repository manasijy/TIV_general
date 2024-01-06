import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import curve_fit

filename = 'equispaced.txt' #'manasijy/TIV_general/equispaced.txt'
with open(filename, 'r') as file:
    # data = file.read()
    data = pd.read_csv(file, sep="\t") 
# sigma = data['column_name']
sigma = data.iloc[:,2]
sigma_array = sigma.to_numpy()
# rhof_expt = (sigma_array/(alpha*G*b))^2
# Ls = 1e-6# %m 
M = 3.01# % 2.96
b = 2.86e-10#  % m
# k1 = 2e8
# k2= 6 #%5-10;
# kL = 150 # % 100-400
alpha = 1/3
G = 26e9 # Pa
sigma_i = 92.35e6 # %Pa
rhof_expt = (sigma_array/(alpha*G*b))^2
# rhof0 = 1e11 #m-2
# Define the system of ODEs
def system(y, t, kL, Ls, k1, k2):
    L, rhof = y
    dL = -kL*(L-Ls)
    drhof = M*((k1/(b*L))- k2*rhof)
    # drhom = (M/b)*(1/Ls-1/L);
    dydt = [dL, drhof]
    return dydt

# Function to fit
def fit_func(t, L0, rhof0, kL, Ls, k1, k2):
    y0 = [L0, rhof0]
    solution = odeint(system, y0, t, args=(kL, Ls, k1, k2,))
    return solution[:, 0]  # replace with the appropriate function of u, v, w

# Time points
t = np.linspace(0, 10, 1000)

# Experimental data
y_exp = np.random.rand(1000)  # replace with your experimental data

# Initial guess for the parameters
# p0 = [1.0, 0.0, 0.0]  # replace with your initial guess
p0 = [100e6, 1e11, 150, 1e-6, 2e8, 6 ]
# Fit the function to the data
popt, pcov = curve_fit(fit_func, t, y_exp, p0)

# Print the optimal parameters
print(popt)
