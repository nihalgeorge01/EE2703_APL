'''
EE2703 Assignment 3
Author: Nihal John George (EE19B131)
'''

from pylab import *
import scipy.special as sp

# Part 4 - g(t;A,B)
def g(t, A, B):
    return A*sp.jn(2,t) + B*t

# Part 2 - loadtxt
data = loadtxt('fitting.dat')

t = data[:,0]                       # time instants
y = data[:,1:]                      # all value columns
N = data.shape[0]                   # time instant count
k = data.shape[1]-1                 # value column count
sigma=logspace(-1,-3,k)             # noise stdev

# Part 3,4 - All data with noise label
A0 = 1.05
B0 = -0.105
j2_vals = g(t, A0, B0)

plot_legend = [r'$\sigma_{' + str(i+1) + r'}=' + '%.3f'%(sigma[i]) + r' $' for i in range(k)]
plot_legend += ['True Value']

figure(0)
plot(t,y, linewidth=1.0)
plot(t,j2_vals, 'k', linewidth=2.0)
xlabel(r'$t \rightarrow$',size=10)
ylabel(r'$f(t) + n \rightarrow$',size=10)
title(r'Q4: Data To Be Fitted To Theory')
legend(plot_legend)
grid(True)
savefig('Q4.png')
show()

# Part 5 - Errorbars
col1 = y[:, 0]

figure(1)
xlabel(r'$t \rightarrow$',size=10)
ylabel(r'$f(t) + n \rightarrow$',size=10)
title(r'Q5: Data With $\sigma = 0.10$ along with Exact Function')
errorbar(t[::5], col1[::5], sigma[0], fmt='ro')
plot(t,j2_vals, 'k', linewidth=2.0)
legend(['f(t)', 'Errorbar'])
grid(True)
savefig('Q5.png')
show()

# Part 6 - Matrix formulation of g(t,A,B) as g(t,1,0)
M = c_[g(t, 1, 0), t]           
p = array([A0,B0])
g_Mp = dot(M,p)
assert array_equal(g_Mp, g(t,A0,B0)), "g(t,A,B) and M.p do not match"

# Part 7 - MSE of data vs g(t,A,B)
A = linspace(0,2,21)
B = linspace(-0.2,0,21)

b = y[:,0]
mse = np.array([[mean(square(b-g(t,i,j))) for j in B] for i in A])

# Part 8 and 9 - Contour plot of MSE and lstsq estimate of A and B
figure(2)
ls_est = lstsq(M,b,rcond=None)
A_est = ls_est[0][0]
B_est = ls_est[0][1]

ind = unravel_index(argmin(mse, axis=None), mse.shape)

plot(A[ind[0]], B[ind[1]], 'ro')
annotate('argmin (%.3f, %.3f)' % (A[ind[0]], B[ind[1]]), (A[ind[0]], B[ind[1]]))

plot(A_est,B_est, 'bo')
annotate('lstsq (%.3f, %.3f)' % (A_est, B_est), (A_est, B_est))

cont = contour(A,B,mse, levels=20)
xlabel(r'A $\rightarrow$',size=10)
ylabel(r'B $\rightarrow$',size=10)
clabel(cont, cont.levels[:5])
title(r'Q8: Contour Plot of $\varepsilon_{ij}$')
grid(True)
savefig('Q8.png')
show()
print("A,B with least MSE:",A_est, B_est)

# Part 10 - Repeat lstsq with each column, plot error in A and B vs noise
temp = [lstsq(M,y[:,i],rcond=None) for i in range(k)]
A = array([i[0][0] for i in temp])
B = array([i[0][1] for i in temp])

A_err = abs(A - A0)
B_err = abs(B - B0)

figure(3)
plot(sigma, A_err, 'ro--', linewidth=0.5, markersize=5.0)
plot(sigma, B_err, 'go--', linewidth=0.5, markersize=5.0)
xlabel(r'Noise Standard Deviation $\rightarrow$',size=10)
ylabel(r'Error $\rightarrow$',size=10)
legend(['A_err', 'B_err'])
title(r'Q10: Variation of Error With Noise')
grid(True)
savefig('Q10.png')
show()

# Part 11 - Replot Part 10 curves using loglog
figure(4)
loglog(sigma, A_err, 'ro--', linewidth=0.5, markersize=5.0)
loglog(sigma, B_err, 'go--', linewidth=0.5, markersize=5.0)
xlabel(r'Noise Standard Deviation $\rightarrow$',size=10)
ylabel(r'Error $\rightarrow$',size=10)
legend(['A_err', 'B_err'])
title(r'Q11: Variation of Error With Noise (loglog)')
grid(True)
savefig('Q11.png')
show()


print("Done")