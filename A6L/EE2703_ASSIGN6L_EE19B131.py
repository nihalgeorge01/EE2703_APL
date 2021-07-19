'''
EE2703 Assignment 6L - Laplace Transform
Author: Nihal John George (EE19B131)
'''

import numpy as np 
import matplotlib.pyplot as plt 
import scipy.signal as sp 

# Q1: Decay = 0.5
freq = 1.5
tc = 2
dec = 1/tc
num = np.poly1d([1,dec])
denom = np.polymul([1,2*dec,dec**2 + freq**2],[1,0,2.25])
H = sp.lti(num,denom)
t,x = sp.impulse(H,None,np.linspace(0,50,501))

plt.figure(1)
plt.plot(t,x)
plt.title(r"Q1: Decaying sinusoidal input (decay = 0.5)")
plt.xlabel(r"t $\rightarrow$")
plt.ylabel(r"x(t) $\rightarrow$")
plt.savefig("Q1_decaying05.png")
#plt.show()

# Q2: Decay = 0.05
freq = 1.5
tc = 20
dec = 1/tc
num = np.poly1d([1,dec])
denom = np.polymul([1,2*dec,dec**2 + freq**2],[1,0,2.25])
H = sp.lti(num,denom)
t,x = sp.impulse(H,None,np.linspace(0,50,501))

plt.figure(2)
plt.plot(t,x)
plt.title(r"Q2: Decaying sinusoidal input (decay = 0.05)")
plt.xlabel(r"t $\rightarrow$")
plt.ylabel(r"x(t) $\rightarrow$")
plt.savefig("Q2_decaying005.png")
plt.show()

# Q3: Simulating using transfer function
H = sp.lti([1],[1,0,2.25])
t = np.linspace(0,50,501)

ct = 3
for freq in np.linspace(1.4,1.6,5):
    f = np.cos(freq*t)*np.exp(-dec*t)
    t,x,svec = sp.lsim(H,f,t) 
    plt.figure(ct)
    plt.plot(t,x)
    plt.title(r"Q3: sp.lsim() simulation (freq = %.3f rad/s, decay = 0.05)" % freq)
    plt.xlabel(r"t $\rightarrow$")
    plt.ylabel(r"x(t) $\rightarrow$")
    plt.savefig("Q3_" + str(ct-2) + "_magphase.png")
    ct += 1
    plt.show()

# Q4: Coupled spring
num = [1,0,2]
denom = [1,0,3,0]
H = sp.lti(num,denom)
t = np.linspace(0,20,201)
t,x = sp.impulse(H, None, t)

num = [2]
denom = [1,0,3,0]
H = sp.lti(num,denom)
t = np.linspace(0,20,201)
t,y = sp.impulse(H, None, t)

plt.figure(ct)
plt.plot(t,x, 'r')
plt.plot(t,y, 'b')
plt.title(r"Q4 : Coupled Spring")
plt.xlabel(r"t $\rightarrow$")
plt.ylabel(r"x(t) $\rightarrow$")
plt.legend(['x','y'])
plt.savefig("Q4_coupled_spring.png")
plt.show()
ct+=1

# Q5: RLC Circuit Transfer Function
H = sp.lti([1],[1e-12, 1e-4, 1])
w,S,phi = H.bode()

plt.figure(ct)
plt.subplot(2,1,1)
plt.semilogx(w,S)
plt.title("Q5: RLC Bode Plot")
plt.xlabel(r"w $\rightarrow$")
plt.ylabel(r"|H| (dB) $\rightarrow$")

plt.subplot(2,1,2)
plt.semilogx(w,phi)
plt.xlabel(r"w $\rightarrow$")
plt.ylabel(r"phase(H) (degrees) $\rightarrow$")
plt.savefig("Q5_magphase.png")
plt.show()
ct+=1

# Q6: Passing input to RLC
val_ct = int(1e5)
t = np.linspace(0,1e-2, val_ct+1)
f = np.cos(1e3 * t) - np.cos(1e6 * t)
t,v,svec = sp.lsim(H,f,t)

plt.figure(ct)
plt.plot(t,v)
plt.title(r"Q6: sp.lsim() simulation of RLC on long timescale (10ms)")
plt.xlabel(r"t $\rightarrow$")
plt.ylabel(r"$v_{o}(t)$ $\rightarrow$")
plt.savefig("Q6_long.png")
ct += 1

plt.figure(ct)
plt.plot(t[:int(30*val_ct/(1e4))],v[:int(30*val_ct/(1e4))])
plt.title(r"Q6: sp.lsim() simulation of RLC on short timescale (30us)")
plt.xlabel(r"t $\rightarrow$")
plt.ylabel(r"$v_{o}(t)$ $\rightarrow$")
plt.savefig("Q6_short.png")
ct += 1
plt.show()