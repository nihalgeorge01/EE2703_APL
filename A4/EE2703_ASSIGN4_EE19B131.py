'''
EE2703 Assignment 4 - Fourier Approximations

Author - Nihal John George (EE19B131)
'''

# Not using pylab since it imports into public space, abs and np.abs collide etc.

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# _e denotes exp(x) functions or variables, _cc denotes cos(cos(x)) functions or variables

# Part 1 - Define functions and plot 
def _e(x):
    # Get exp(x)
    return np.exp(x)

def _e_prd(x):
    # Get exp(x) over 0 -> 2 PI, periodically extended using modulus
    return np.exp(x%(2*np.pi))

def _cc(x):
    # Get cos(cos(x))
    return np.cos(np.cos(x))

def _u(x, k, f):
    # Get cosine term
    return f(x)*np.cos(k*x)

def _v(x, k, f):
    # Get sine term
    return f(x)*np.sin(k*x)

def get_A_mat(x,k):
    # Matrix of functions for least squares
    T = x.shape[0]
    A = np.zeros((T,k))

    A[:,0] = 1
    for i in range(1,k//2+1):
        A[:,2*i - 1] = np.cos(i*x)
        A[:,2*i] = np.sin(i*x)

    return A

def get_lstsq(x, k, f):
    # Applies function f on to get f(x) = b. Then gets A matrix using get_A_mat(x,k)
    # Finally solves Ax=b using least squares method

    b = f(x)
    A = get_A_mat(x,k)
    
    c = np.linalg.lstsq(A,b, rcond=None)[0]
    return c

FUNC_MAP = {'exp': _e, 'coscos': _cc}

# Part 2 - First n coeffs function
def eval_coeffs(n, funcstr):
    # Evaluate first n Fourier coefficients (n//2 cosine, n//2 sine, 1 DC) of
    # function denoted by func_str

    func = FUNC_MAP[funcstr]

    a0 = quad(_u, 0, 2*np.pi, args=(0, func))[0]/(2*np.pi)
    a, b = np.empty(n//2), np.empty(n//2)

    for i in range(n//2):
        a[i] = quad(_u, 0, 2*np.pi, args=(i+1, func))[0]/np.pi
        b[i] = quad(_v, 0, 2*np.pi, args=(i+1, func))[0]/np.pi
    
    mixed = np.empty(a.size + b.size + 1, dtype=a.dtype)
    mixed[0] = a0
    mixed[1::2], mixed[2::2] = a, b             # interleave a and b values

    return mixed

x = np.linspace(-2*np.pi, 4*np.pi, 400)
e_vec = _e(x)
e_prd_vec = _e_prd(x)
cc_vec = _cc(x)

plt.figure(1)
plt.semilogy(x,e_vec,'r-')
plt.semilogy(x, e_prd_vec, 'b-')
plt.xlabel(r'$x \rightarrow$',size=10)
plt.ylabel(r'$exp(x) \rightarrow$',size=10)
plt.title(r'exp(x)')
plt.legend(["Actual", "Periodic Extension"])
plt.grid(True)
plt.savefig('Q1_exp.png')
plt.show()

plt.figure(2)
plt.plot(x, cc_vec, 'r-')
plt.xlabel(r'$x \rightarrow$',size=10)
plt.ylabel(r'$cos(cos(x)) \rightarrow$',size=10)
plt.title(r'cos(cos(x))')
plt.grid(True)
plt.savefig('Q1_cos_cos.png')
plt.show()

# Part 2 - First 51 coeffs

mixed_e = eval_coeffs(51, 'exp')
mixed_cc = eval_coeffs(51, 'coscos')
x = range(51)

# Part 3 - Plotting the coeffs

plt.figure(3)
plt.semilogy(x, np.abs(mixed_e), 'ro')
plt.xlabel(r'Coefficient index $\rightarrow$',size=10)
plt.ylabel(r'Coefficient value $\rightarrow$',size=10)
plt.title(r'exp(x) Coefficients (Semilog)')
plt.grid(True)
plt.savefig('Q3_exp_semilog.png')
plt.show()

plt.figure(4)
plt.loglog(x, np.abs(mixed_e), 'ro')
plt.xlabel(r'Coefficient index $\rightarrow$',size=10)
plt.ylabel(r'Coefficient value $\rightarrow$',size=10)
plt.title(r'exp(x) Coefficients (loglog)')
plt.grid(True)
plt.savefig('Q3_exp_loglog.png')
plt.show()

plt.figure(5)
plt.semilogy(x, np.abs(mixed_cc), 'ro')
plt.xlabel(r'Coefficient index $\rightarrow$',size=10)
plt.ylabel(r'Coefficient value $\rightarrow$',size=10)
plt.title(r'cos(cos(x)) Coefficients (Semilog)')
plt.grid(True)
plt.savefig('Q3_cos_cos_semilog.png')
plt.show()

plt.figure(6)
plt.loglog(x, np.abs(mixed_cc), 'ro')
plt.xlabel(r'Coefficient index $\rightarrow$',size=10)
plt.ylabel(r'Coefficient value $\rightarrow$',size=10)
plt.title(r'cos(cos(x)) Coefficients (loglog)')
plt.grid(True)
plt.savefig('Q3_cos_cos_loglog.png')
plt.show()

# Part 4 - Least Squares Approach
x = np.linspace(0, 2*np.pi, 401)[:-1]
c_e = get_lstsq(x, 51, _e)
c_cc = get_lstsq(x, 51, _cc)
ind = np.array(range(51))               # coefficient index vector

# Part 5 - Plot Least Squares coefficients
plt.figure(7)
plt.semilogy(ind, np.abs(mixed_e), 'ro')
plt.semilogy(ind, np.abs(c_e), 'go')
plt.xlabel(r'Coefficient index $\rightarrow$',size=10)
plt.ylabel(r'Coefficient value $\rightarrow$',size=10)
plt.title(r'exp(x) Coefficients (Semilog)')
plt.legend(['Integral', 'Least Squares'])
plt.grid(True)
plt.savefig('Q5_exp_semilog.png')
plt.show()

plt.figure(8)
plt.loglog(ind, np.abs(mixed_e), 'ro')
plt.loglog(ind, np.abs(c_e), 'go')
plt.xlabel(r'Coefficient index $\rightarrow$',size=10)
plt.ylabel(r'Coefficient value $\rightarrow$',size=10)
plt.title(r'exp(x) Coefficients (Loglog)')
plt.legend(['Integral', 'Least Squares'])
plt.grid(True)
plt.savefig('Q5_exp_loglog.png')
plt.show()

plt.figure(9)
plt.semilogy(ind, np.abs(mixed_cc), 'ro')
plt.semilogy(ind, np.abs(c_cc), 'go')
plt.xlabel(r'Coefficient index $\rightarrow$',size=10)
plt.ylabel(r'Coefficient value $\rightarrow$',size=10)
plt.title(r'cos(cos(x)) Coefficients (Semilog)')
plt.legend(['Integral', 'Least Squares'])
plt.grid(True)
plt.savefig('Q5_cos_cos_semilog.png')
plt.show()

plt.figure(10)
plt.loglog(ind, np.abs(mixed_cc), 'ro')
plt.loglog(ind, np.abs(c_cc), 'go')
plt.xlabel(r'Coefficient index $\rightarrow$',size=10)
plt.ylabel(r'Coefficient value $\rightarrow$',size=10)
plt.title(r'cos(cos(x)) Coefficients (Loglog)')
plt.legend(['Integral', 'Least Squares'])
plt.grid(True)
plt.savefig('Q5_cos_cos_loglog.png')
plt.show()

# Part 6 - Get deviation b/w least squares and integration
maxerr_e = np.max(np.abs(mixed_e-c_e))
maxerr_cc = np.max(np.abs(mixed_cc-c_cc))

print("Deviation in exp(x):", maxerr_e)
print("Deviation in cos(cos(x)):", maxerr_cc)

# Part 7 - Plot function estimate by least squares and compare with true function
x = np.linspace(0, 2*np.pi, 400)

est_e = np.matmul(get_A_mat(x, 51),c_e)             # least square estimates of the functions
est_cc = np.matmul(get_A_mat(x, 51),c_cc)
e_prd_vec = _e_prd(x)                               # actual values of the functions
cc_vec = _cc(x)

plt.figure(1)
plt.semilogy(x,est_e,'go')
plt.semilogy(x, e_prd_vec, 'r-')
plt.xlabel(r'$x \rightarrow$',size=10)
plt.ylabel(r'$exp(x) \rightarrow$',size=10)
plt.title(r'exp(x) Estimate')
plt.legend(["Least Square", "Actual", "Periodic Extension"])
plt.grid(True)
plt.savefig('Q7_exp.png')
plt.show()

plt.figure(2)
plt.plot(x, est_cc, 'go')
plt.plot(x, cc_vec, 'r-')
plt.xlabel(r'$x \rightarrow$',size=10)
plt.ylabel(r'$cos(cos(x)) \rightarrow$',size=10)
plt.title(r'cos(cos(x)) Estimate')
plt.legend(["Least Square", "Actual"])
plt.grid(True)
plt.savefig('Q7_cos_cos.png')
plt.show()
