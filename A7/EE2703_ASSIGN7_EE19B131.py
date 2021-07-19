'''
EE2703 Assignment 7 - Analysis of Circuits Using Sympy

Author: Nihal John George (EE19B131)
'''

import sympy as sm 
import matplotlib.pyplot as plt 
import numpy as np
import scipy.signal as sp 

PI = np.pi
fig_ct = 1

def lowpass(R1,R2,C1,C2,G,Vi):
    A = sm.Matrix([
                    [0,0,1,-1/G] \
                    , [-1/(1+s*R2*C2),1,0,0] \
                    , [0,-G,G,1] \
                    , [-(1/R1)-(1/R2)-s*C1,1/R2,0,s*C1]
                 ])
    
    b = sm.Matrix([0,0,0,-Vi/R1])
    V = A.inv()*b
    return (A,b,V)

def highpass(R1,R3,C1,C2,G,Vi):
    A = sm.Matrix([
                    [0,0,1,-1/G] \
                    , [-(s*C2*R3)/(1+s*C2*R3),1,0,0] \
                    , [0,-G,G,1] \
                    , [-(1/R1)-s*C2-s*C1,s*C2,0,1/R1]
                 ])
    
    b = sm.Matrix([0,0,0,-Vi*s*C1])
    V = A.inv()*b
    return (A,b,V)

def sm_to_sp(expr):
    # Converts Sympy expression to Scipy.signal compatible numerator and denominator

    num, denom = sm.simplify(expr).as_numer_denom()
    num, denom = sm.poly(num).all_coeffs(), sm.poly(denom).all_coeffs()
    num = [float(i) for i in num]
    denom = [float(i) for i in denom]
    return sp.lti(num,denom)
    

x,y,z = sm.symbols("x y speed")
s = sm.symbols('s')

# Q1: Step Response
A,b,V = lowpass(10000,10000,1e-9,1e-9,1.586,1/s)
Vo = V[3]
H = sm_to_sp(Vo)
t = np.linspace(0,2e-4,int(1e5)+1)
t, out = sp.impulse(H,None,t)

plt.figure(fig_ct)
plt.plot(t,out)
plt.title("Q1: Step Response")
plt.grid(True)
plt.xlabel(r'$t \rightarrow$')
plt.ylabel(r'$V_{o}(t) \rightarrow$')
plt.savefig(str(fig_ct)+'.png')
plt.show()
fig_ct+=1

# Q2: Sinusoidal inputs
A,b,V = highpass(10000,10000,1e-9,1e-9,1.586,1)
H = sm_to_sp(V[3])
t = np.linspace(0,0.00001,int(1e5)+1)
vi = np.sin(2e3 * PI * t) + np.cos(2e6 * PI * t)
t,out,svec = sp.lsim(H,vi,t)

plt.figure(fig_ct)
plt.plot(t,out)
plt.title("Q2: Sinusoidal Input Response")
plt.grid(True)
plt.xlabel(r'$t \rightarrow$')
plt.ylabel(r'$V_{o}(t) \rightarrow$')
plt.savefig(str(fig_ct)+'.png')
plt.show()
fig_ct+=1

# Q3: New circuit transfer function plot
A,b,V = highpass(10000,10000,1e-9,1e-9,1.586,1)
Vo = V[3]
w = np.logspace(0,8,801)
ss = 1j*w
hf = sm.lambdify(s,Vo,'numpy')
v = hf(ss)

plt.figure(fig_ct)
plt.loglog(w,abs(v),lw=2)
plt.title("Q3: Highpass Bode Magnitude Plot")
plt.grid(True)
plt.xlabel(r'$\omega \rightarrow$')
plt.ylabel(r'$|H(j\omega)| \rightarrow$')
plt.savefig(str(fig_ct)+'.png')
plt.show()
fig_ct+=1

# Q4: Response to damped sinusoid
H = sm_to_sp(Vo)
t_low = np.linspace(0,1,int(1e5)+1)
t_high = np.linspace(0,1e-5,int(1e5)+1)
vi_high = np.sin(1e7 * PI * t_high)*np.exp(-1.5e5*t_high)
vi_low = np.sin(1e2 * PI * t_low)*np.exp(-1.5*t_low)
t_high,out_high,svec = sp.lsim(H,vi_high,t_high)
t_low,out_low,svec = sp.lsim(H,vi_low,t_low)

plt.figure(fig_ct)
plt.plot(t_low,out_low)
plt.title("Q4: Highpass Decaying Low-Freq-Sinusoidal Input Response")
plt.grid(True)
plt.xlabel(r'$t \rightarrow$')
plt.ylabel(r'$V_{o}(t) \rightarrow$')
plt.savefig(str(fig_ct)+'.png')
plt.show()
fig_ct+=1

plt.figure(fig_ct)
plt.plot(t_high,out_high)
plt.title("Q4: Highpass Decaying High-Freq-Sinusoidal Input Response")
plt.grid(True)
plt.xlabel(r'$t \rightarrow$')
plt.ylabel(r'$V_{o}(t) \rightarrow$')
plt.savefig(str(fig_ct)+'.png')
plt.show()
fig_ct+=1

# Q5: Step Response
A,b,V = highpass(10000,10000,1e-9,1e-9,1.586,1/s)
Vo = V[3]
H_step = sm_to_sp(Vo)
t = np.linspace(0,0.001,int(1e5)+1)
t, out_step = sp.impulse(H_step,None,t)

plt.figure(fig_ct)
plt.plot(t,out_step)
plt.title("Q5: Highpass Step Response")
plt.grid(True)
plt.xlabel(r'$t \rightarrow$')
plt.ylabel(r'$V_{o}(t) \rightarrow$')
plt.savefig(str(fig_ct)+'.png')
plt.show()
fig_ct+=1