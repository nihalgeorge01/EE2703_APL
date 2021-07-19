'''
Assignment 5 - The Resistor Problem
Author - Nihal John George (EE19B131)
'''

# Part 1 - Imports
from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
import seaborn as sns
from argparse import ArgumentParser

def leastsq_fit(x,y):
    
    # reshaping to make it concat compatible
    x = x.reshape((-1,1))
    one_arr = np.ones((x.shape[0],1))

    # left side matrix of 2 columns : 1st column with 1's, 
    # 2nd column with coeffs
    xvar = concatenate((one_arr,x), axis=1)
    
    logA, B = lstsq(xvar, log(y), rcond=None)[0]
    return exp(logA), B

# Part 2 - Parse Arguments
ap = ArgumentParser()

ap.add_argument('--Nx', default=25, type=int, help='X Resolution')
ap.add_argument('--Ny', default=25, type=int, help='Y Resolution')
ap.add_argument('--radius', default=8, type=int,help='Radius in indices')
ap.add_argument('--Niter', default=1500, type=int, help='Number of Iterations')

argv = ap.parse_args()

Nx, Ny, Nr, Niter = argv.Nx, argv.Ny, argv.radius, argv.Niter

# Part 3 - Initialize Potential Array
phi = zeros((Ny,Nx), dtype='float64')

x = linspace(0,Nx-1,Nx) - Nx//2             # Subtraction with Nx//2 to centre the distribution
y = linspace(0,Ny-1,Ny) - Ny//2             
Y, X = meshgrid(y,x)
ii = np.where(X*X + Y*Y <= Nr*Nr)
phi[ii] = 1.0                               # set aspect ratio to 1
phi_orig = phi.copy()
ii_T = np.array(ii).T

# Contour plot with 1V nodes marked red
fig1 = figure(1)
contourf(Y,X,phi)
plot(ii_T[:,0]-Ny//2, ii_T[:,1]-Nx//2, 'ro')            # Subtraction with Nx and Ny to centre the plot to (0,0)
ax = fig1.axes
ax[0].set_aspect(1)
title("Diagram of Plate and Electrode")
xlabel(r"x $\rightarrow$")
ylabel(r"y $\rightarrow$")
savefig('electrode.png')
show()

# Part 4,5,6 - Perform Iteration, Update Potential 
# and Enforce Boundary Conditions
oldphi = phi.copy()
errors = np.zeros(Niter)

for k in range(Niter):
    oldphi = phi.copy()
    phi[1:-1, 1:-1] = 0.25*(oldphi[1:-1, 0:-2] + \
                            oldphi[1:-1, 2:] + \
                            oldphi[0:-2, 1:-1] + \
                            oldphi[2:, 1:-1] \
                            )
    
    # Boundary update
    phi[1:-1, 0] = phi[1:-1,1]          # Left Boundary
    phi[1:-1, -1] = phi[1:-1, -2]       # Right Boundary
    phi[0, 1:-1] = phi[1, 1:-1]         # Top Boundary
    phi[-1, 1:-1] = 0                   # Bottom Boundary (Ground)

    phi[ii] = 1.0                       # Wire Boundary

    errors[k] = (abs(phi-oldphi)).max()

# Part 7 - Fitting with straight lines
ind = array(range(0,Niter))
ind1 = array(range(0,Niter))
ind2 = ind1[500:]
A1,B1 = leastsq_fit(ind1, errors)               # fit with all points
A2,B2 = leastsq_fit(ind2, errors[500:])         # fit with all except first 500

fit1 = A1*exp(B1*ind1)
fit2 = A2*exp(B2*ind2)

# semilog plots
figure(2)
semilogy(ind[::50],errors[::50], 'ro--', markersize=12)
semilogy(ind1[::50], fit1[::50], 'go--', markersize=8)
semilogy(ind2[::50], fit2[::50], 'bo--', markersize=4)
legend(('error', 'fit all', 'fit 500+'))
title("Error and Error-Fitting Plots")
xlabel(r"Iteration Count $\rightarrow$")
ylabel(r"Error $\rightarrow$")
savefig('errors.png')
show()

# loglog plot
figure(3)
loglog(ind[::50],errors[::50], 'ro--', markersize=12)
title("Error Plot (LogLog)")
xlabel(r"Iteration Count $\rightarrow$")
ylabel(r"Error $\rightarrow$")
savefig('error_loglog.png')
show()

# Part 8 - Error estimation
maxerr_1 = abs(-A1*exp(B1*(Niter+0.5)))
maxerr_2 = abs(-A2*exp(B2*(Niter+0.5)))

print("Error Bound for fit 1:", maxerr_1)
print("Error Bound for fit 2:", maxerr_2)

# Part 9 - Surface Plot
fig3 = figure(3)
ax = p3.Axes3D(fig3)
surf = ax.plot_surface(X,Y, phi.T[:, ::-1], rstride=1, cstride=1, cmap=cm.jet)      # Reverse the columns of phi.T to get correct coordinates on the plot
xlabel(r"X $\rightarrow$")
ylabel(r"Y $\rightarrow$")
ax.set_zlabel(r"Potential $\rightarrow$")
ax.set_title("Surface Plot of Potential")
savefig('surface.png')
show()

# Part 10 - Contour Plot
fig4 = figure(4)
contourf(Y,X[::-1], phi)
plot(ii_T[:,0]-Ny//2, ii_T[:,1]-Nx//2, 'ro')                # Subtraction with Nx and Ny to centre the plot to (0,0)
ax = fig4.axes
ax[0].set_aspect(1)
title("Contour Plot of Potential")
xlabel(r"x $\rightarrow$")
ylabel(r"y $\rightarrow$")
savefig('contour_last.png')
show()

# Part 11 - Vector Plot of Current 
Jx = np.zeros(phi.shape)
Jy = np.zeros(phi.shape)
Jx[:,1:-1] = (phi[:,0:-2]-phi[:,2:])/2
Jy[1:-1,:] = (phi[2:, :]-phi[0:-2,:])/2

J_tot_sq = Jx**2 + Jy**2

fig5 = figure(5)
quiver(Y,X[::-1],Jx,Jy, scale=6)
plot(ii_T[:,0]-Ny//2, ii_T[:,1]-Nx//2, 'ro')    # Subtraction with Nx and Ny to centre the plot to (0,0)
ax = fig5.axes
ax[0].set_aspect(1)                             # Set aspect ratio to 1
title("Vector Plot of Current")
xlabel(r"x $\rightarrow$")
ylabel(r"y $\rightarrow$")
savefig('quiver.png')
show()

# Part 12 (Optional) - Heat Plot
# Taking heat and electrical conductivity as 1
# del^2(T) = J^2

temp = np.zeros(phi.shape) + 300.0

oldtemp = temp.copy()

for k in range(Niter):
    oldtemp = temp.copy()
    temp[1:-1, 1:-1] = 0.25*(  oldtemp[1:-1, 0:-2] \
                            + oldtemp[1:-1, 2:] \
                            + oldtemp[0:-2, 1:-1] \
                            + oldtemp[2:, 1:-1] \
                            + 10*J_tot_sq[1:-1,1:-1] \
                            )
    
    # Boundary update
    temp[1:-1, 0] = temp[1:-1,1]
    temp[1:-1, -1] = temp[1:-1, -2]
    temp[0, 1:-1] = temp[1, 1:-1]
    temp[-1, 1:-1] = 300.0 # connected to ground

    temp[ii] = 300.0

figure(7)
ax = sns.heatmap(temp, xticklabels=x, yticklabels=y[::-1], cmap='coolwarm')
ax.set_title('Heat Map')
xlabel(r"x $\rightarrow$")
ylabel(r"y $\rightarrow$")
savefig('temp.png')
show()