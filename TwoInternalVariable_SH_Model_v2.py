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
sigma_array_dis = sigma_array-sigma_array[0]+M*alpha*G*b #*np.sqrt(1e12)
rhof_expt = (sigma_array_dis*1e6/(M*alpha*G*b))**2

# rhof0 = 1e11 #m-2
# Define the system of ODEs
def system(y, t, kL, Ls, k1, k2):
    L, rhof, rhom = y
    dL = -kL*(L-Ls)
    drhof = M*((k1/(b*L))- k2*rhof)
    drhom = (M/b)*(1/Ls-1/L)
    dydt = [dL, drhof, drhom]
    return dydt

# Function to fit

def fit_func(t, L0, rhof0, rhom0, kL, Ls, k1, k2):
    y0 = [L0, rhof0, rhom0]
    solution = odeint(system, y0, t, args=(kL, Ls, k1, k2))
    return solution[:, 1]

# Time points
# t = np.linspace(0, 10, 1000)
t = strain

# Initial guess for the parameters
p0 = [40e-6, 1e12, 1e12, 1, 1e-6, 1e8, 6 ]
# Fit the function to the data
L_U_bounds = ([39e-6, 9e11, 9e11, 0.001, 1e-7, 1e6, 1],
           [40e-6, 1.1e12, 1.1e12, 100, 30e-6, 1e10, 200])
popt, pcov = curve_fit(fit_func, t, rhof_expt, p0, bounds=L_U_bounds)
# Print the optimal parameters
print('L0=',"{:.6e}".format(popt[0])),print('rhof0=',"{:.6e}".format(popt[1]))
print('rhom0=',"{:.6e}".format(popt[2]))
print('kL=',"{:.6e}".format(popt[3])), print('Ls=',"{:.6e}".format(popt[4]))
print('k1=',"{:.6e}".format(popt[5])),print('k2=',"{:.6e}".format(popt[6]))
## comparing the result against experimental data
import matplotlib.pyplot as plt

# Time points
# t = np.linspace(0, 10, 1000)

# Optimal parameters from curve_fit
L0_opt, rhof0_opt, rhom0_opt, kL_opt, Ls_opt, k1_opt, k2_opt = popt
rhof_fit = fit_func(t, L0_opt, rhof0_opt, rhom0_opt, kL_opt, Ls_opt, k1_opt, k2_opt)

# Plot the experimental data
plt.scatter(t, rhof_expt, label='Experimental data')

# Plot the fitted curve
plt.plot(t, rhof_fit, label='Fitted curve', color='red')

# Add a legend
plt.legend()

# Show the plot
plt.show()


# Initial conditions with optimal parameters
y0_opt = [L0_opt, rhof0_opt, rhom0_opt]

# Solve the ODEs with the optimal parameters
solution_opt = odeint(system, y0_opt, t, args=(kL_opt, Ls_opt, k1_opt, k2_opt))


# Extract the solutions for u, v, and w
L_opt = solution_opt[:, 0]
rhof_opt = solution_opt[:, 1]
rhom_opt = solution_opt[:, 2]
plt.plot(t, rhof_opt, label='L', color='red')
plt.show()
plt.plot(t, L_opt, label='rho_f', color='green')
plt.plot(t, rhom_opt, label='rho_m', color='blue')
plt.show()


# # Define your fit function with only the parameters to optimize
# def fit_func(t, param1, param2, param3):
#     u0, v0, w0 = 40e-6, 1e10, 1e10  # fixed parameters
#     a, b, c = param1, param2, param3  # parameters to optimize
#     y0 = [u0, v0, w0]
#     solution = odeint(system, y0, t, args=(a, b, c))
#     return solution[:, 0]  # replace with the appropriate function of u, v, w

# # Initial guess for the parameters to optimize
# p0 = [0.001, 1e-7, 1e6]  # replace with your initial guess

# # Fit the function to the data
# popt, pcov = curve_fit(fit_func, t, y_exp, p0)
