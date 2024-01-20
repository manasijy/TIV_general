import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import curve_fit

filename = 'C:/Users/MKY/Desktop/equispaced.txt'
 #'manasijy/TIV_general/equispaced.txt'
with open(filename, 'r') as file:
    # data = file.read()
    data = pd.read_csv(file, sep="\t") 
# sigma = data['column_name']
strain = data.iloc[:,0]
sigma = data.iloc[:,1]
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
# sigma_i = 92.35e6 # %Pa
rhof_expt = (sigma_array*1e6/(M*alpha*G*b))**2
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
    return solution[:, 1]

# Time points
# t = np.linspace(0, 10, 1000)
t = strain

# Initial guess for the parameters
p0 = [40e-6, 1e12, 1, 1e-6, 1e8, 6 ]
# Fit the function to the data
popt, pcov = curve_fit(fit_func, t, rhof_expt, p0)
# Print the optimal parameters
print('L0=',"{:.6e}".format(popt[0])),print('rhof0=',"{:.6e}".format(popt[1]))
print('kL=',"{:.6e}".format(popt[2])), print('Ls=',"{:.6e}".format(popt[3]))
print('k1=',"{:.6e}".format(popt[4])),print('k2=',"{:.6e}".format(popt[5]))
## comparing the result against experimental data
import matplotlib.pyplot as plt

# Time points
# t = np.linspace(0, 10, 1000)

# Use the optimal parameters to generate the fitted curve
L0, rhof0, kL, Ls, k1, k2 = popt
rhof_fit = fit_func(t, L0, rhof0, kL, Ls, k1, k2)

# Plot the experimental data
plt.scatter(t, rhof_expt, label='Experimental data')

# Plot the fitted curve
plt.plot(t, rhof_fit, label='Fitted curve', color='red')

# Add a legend
plt.legend()

# Show the plot
plt.show()
