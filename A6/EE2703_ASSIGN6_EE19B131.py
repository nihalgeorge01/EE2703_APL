'''
Assignment 6: Tubelight Simulation
Takes about 5 seconds to run with default params

Author: Nihal John George (EE19B131)

Run with defaults
python EE2703_ASSIGN6_EE19B131.py

Run with custom values
python EE2703_ASSIGN6_EE19B131.py --n 100 --M 5 --nk 500 --u0 7 --p 0.5 --Msig 0.1
'''

import sys
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib as mpl
from argparse import ArgumentParser

# This line prevents error due to too many values to plot in phase space
mpl.rcParams['agg.path.chunksize'] = 10000

# Using argparse instead of sys.argv
# since it can handle argcount and type mismatch better
ap = ArgumentParser()

ap.add_argument('--n', default=100, type=int, help='Spatial Grid Size')
ap.add_argument('--M', default=5, type=int, help='Number of electrons injected per turn')
ap.add_argument('--nk', default=500, type=int,help='Number of turns to simulate')
ap.add_argument('--u0', default=7.0, type=float, help='Threshold Velocity')
ap.add_argument('--p', default=0.5, type=float, help='Probability that ionization will occur')
ap.add_argument('--Msig', default=0.1, type=float, help='Standard deviation of electron injection distribution')

argv = ap.parse_args()
n, M, nk, u0, p, Msig = argv.n, argv.M, argv.nk, argv.u0, argv.p, argv.Msig

# Initialise position displacement and velocity holders
xx = np.zeros(n*M)
u = np.zeros(n*M)
dx = np.zeros(n*M)

# Initialise intensity, position and velocity accumulators
I = []
X = []
V = []

# Iterate over timesteps
# Initialise ii here and once per iteration
ii = np.where(xx>0)
for k in range(1,nk+1):
    # Update kinematic variables due to field acceleration
    dx[ii] = u[ii] + 0.5
    xx[ii] += dx[ii]
    u[ii] += 1

    # Reset electrons which hit anode to free slots
    hit_el = np.where(xx>=n)
    xx[hit_el] = 0
    dx[hit_el] = 0
    u[hit_el] = 0

    # Ionize electrons with high energy with probability p
    high_el = np.array(np.where(u>=u0)[0])   # kk
    ion_ind = np.array(np.where(np.random.rand(len(high_el))<=p)[0])     #ll
    ion_el = high_el[ion_ind]   # kl
    
    # Set ionized electrons to a probable collision location, and rest
    u[ion_el] = 0
    rho = np.random.rand(len(dx[ion_el]))
    # xx[ion_el] -= dx[ion_el]*rho
    # u[ion_el] -= rho
    
    # Formula taking into account s=ut+0.5at^2 given below
    xx[ion_el] -= dx[ion_el]*(1-np.divide((u[ion_el]*rho) + 0.5*(rho**2), u[ion_el]+0.5))
    u[ion_el] -= (1-rho)

    # Get number of photons at each x (which is also number of electrons at each x)
    I.extend(xx[ion_el].tolist())

    # Inject m new electrons
    m = max(int(np.random.randn()*Msig + M),0)

    # Assign injected electrons to free space, taking care of limited space
    free = np.where(xx==0)
    free = free[:m]
    xx[free] = 1

    # Find existing electrons and note down position and velocity in accumulators
    ii = np.where(xx>0)[0]
    X.extend(xx[ii].tolist())
    V.extend(u[ii].tolist())

# Position Histogram
plt.figure(0)
bin_edges = [i for i in range(n+1)]
ct, bins, _ = plt.hist(X, bins=bin_edges)
plt.title("Position Histogram")
plt.xlabel(r"x $\rightarrow$")
plt.ylabel(r"Count $\rightarrow$")
plt.savefig("pos_hist__" + \
            "_n_" + ''.join(str(n).split('.')) + \
            "_M_" + ''.join(str(M).split('.')) + \
            "_nk_" + ''.join(str(nk).split('.')) + \
            "_u0_" + ''.join(str(u0).split('.')) + \
            "_p_" + ''.join(str(p).split('.')) + \
            "_Msig_" + ''.join(str(Msig).split('.')) + \
            "_.png")

# Intensity Histogram
plt.figure(1)
plt.hist(I, bins=bin_edges)
plt.title("Intensity Histogram")
plt.xlabel(r"x $\rightarrow$")
plt.ylabel(r"Intensity $\rightarrow$")
plt.savefig("int_hist__" + \
            "_n_" + ''.join(str(n).split('.')) + \
            "_M_" + ''.join(str(M).split('.')) + \
            "_nk_" + ''.join(str(nk).split('.')) + \
            "_u0_" + ''.join(str(u0).split('.')) + \
            "_p_" + ''.join(str(p).split('.')) + \
            "_Msig_" + ''.join(str(Msig).split('.')) + \
            "_.png")

# Phase Space (X vs V)
plt.figure(2)
plt.scatter(X,V)
plt.title("Phase Space")
plt.xlabel(r"X $\rightarrow$")
plt.ylabel(r"V $\rightarrow$")
plt.savefig("phase__" + \
            "_n_" + ''.join(str(n).split('.')) + \
            "_M_" + ''.join(str(M).split('.')) + \
            "_nk_" + ''.join(str(nk).split('.')) + \
            "_u0_" + ''.join(str(u0).split('.')) + \
            "_p_" + ''.join(str(p).split('.')) + \
            "_Msig_" + ''.join(str(Msig).split('.')) + \
            "_.png")

# Print table of intensity bins and counts
xpos = 0.5*(bins[0:-1]+bins[1:])
print("Intensity Data")
print("xpos \t count")

for i in range(len(ct)):
    print(xpos[i], '\t', ct[i])

# Show all figures at once
plt.show()
