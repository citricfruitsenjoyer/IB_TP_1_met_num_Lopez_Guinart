import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
plt.style.use('ggplot')
# from plotly import express as px
# from scipy import signal as sn
# from scipy.constants import Rydberg as R
from scipy import optimize as o
plt.rcParams.update({
    "figure.figsize": (12, 9),
    "axes.titlesize": 24,
    "axes.labelsize": 18,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 15})


# T = pd.read_parquet('datos/cubo.parquet')
T = pd.read_csv('datos/cubo.csv',names=['T'])
# T.to_parquet('datos/cubo.parquet', compression='gzip')

T = 250*T/np.max(T)

n = int(round(len(T)**(1/3)));n

d = 1/(n+1)

mid = int(1/(2*d))   # ≈ 128

i = j = k = mid

idx_i = lambda i: i + n*(j-1) + n**2*(k-1)
idx_j = lambda j: i + n*(j-1) + n**2*(k-1)
idx_k = lambda k: i + n*(j-1) + n**2*(k-1)

T_x = np.ones(n)
T_y = np.ones(n)
T_z = np.ones(n)

# (x, 1/2, 1/2)
for ii in range(n):
    T_x[ii] = T['T'].iloc[idx_i(ii+1) - 1]

# (1/2, y, 1/2)
for jj in range(n):
    T_y[jj] = T['T'].iloc[idx_j(jj+1) - 1]

# (1/2, 1/2, z)
for kk in range(n):
    T_z[kk] = T['T'].iloc[idx_k(kk+1) - 1]
    
    
# reconstruir el cubo desde el vector
T_vals = T['T'].values
T_cube = T_vals.reshape((n, n, n), order='F')  # MUY importante

# coordenadas físicas
x = np.linspace(d, 1-d, n)
y = np.linspace(d, 1-d, n)
z = np.linspace(d, 1-d, n)

X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# flatten para scatter
Xf = X.flatten()
Yf = Y.flatten()
Zf = Z.flatten()
Tf = T_cube.flatten()

# plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

p = ax.scatter(Xf, Yf, Zf, c=Tf)

fig.colorbar(p)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('Distribución de temperatura en el cubo')

plt.show()