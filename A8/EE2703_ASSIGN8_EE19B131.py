from pylab import *

# Function for dft since it is used ~7 times
def dft(func, ti, tf, N, limx, ylabel1, ylabel2, title_str, path, thresh=0, graph=True):
    t = linspace(ti,tf,N+1)[:-1]
    y = func(t)

    # Shift negative time half to front and then do fft, then shift back
    # This does not matter for periodic even or periodic odd functions
    Y = fftshift(fft(ifftshift(y)))/N

    # Scale the frequenices to show sensible values
    w = (float(N)/(tf-ti))*(linspace(-pi, pi, N+1)[:-1])

    if graph:
        figure()
        
        subplot(2,1,1)
        plot(w,abs(Y),lw=2)
        xlim([-limx,limx])
        ylabel(ylabel1,size=16)
        title(title_str)
        grid(True)

        subplot(2,1,2)
        
        ii = where(abs(Y)>=thresh)
        plot(w,angle(Y),'ro', ms=2, lw=2)
        plot(w[ii],angle(Y[ii]),'go',ms=10, lw=5)
        xlim([-limx,limx])
        ylabel(ylabel2,size=16)
        xlabel(r"$k$",size=16)
        legend(['too small (<1e-3)', 'valid (>1e-3)'])
        grid(True)
        
        savefig(path)
        show()
    
    return Y

#'''
# FFT and IFFT inverse relation check
x = rand(100)
X = fft(x)
y = ifft(X)
print(c_[x,y])
print("Max error after FFT+IFFT:", abs(x-y).max())

# Sine (incorrect)
x = linspace(0,2*pi,128)
y = sin(5*x)
Y = fft(y)

subplot(2,1,1)
plot(abs(Y), lw=2)
ylabel(r"$|Y|$", size=16)
title(r"Spectrum of $\sin(5t)$")
grid(True)

subplot(2,1,2)
plot(unwrap(angle(Y)), lw=2)
ylabel(r"$\angle Y$", size=16)
xlabel(r"$k$", size=16)
grid(True)
savefig("sine_1.png")
show()

# Sine (correct)
func = lambda x : sin(5*x)
dft(func, 0, 2*pi, 128, 10, \
    r"$|Y|$", r"$\angle Y$", r"Spectrum of $\sin(5t)$", \
    "sine_2.png", 1e-3)

# AM (low res)
func = lambda x : (1+0.1*cos(x))*cos(10*x)
dft(func, 0, 2*pi, 128, 15, \
    r"$|Y|$", r"$\angle Y$", r"Spectrum of $\left(1+0.1\cos\left(t\right)\right)\cos\left(10t\right)$", \
    "AM_lo.png", 1e-3)


# AM (hi res)
func = lambda x : (1+0.1*cos(x))*cos(10*x)
dft(func, -4*pi, 4*pi, 512, 15, \
    r"$|Y|$", r"$\angle Y$", r"Spectrum of $\left(1+0.1\cos\left(t\right)\right)\cos\left(10t\right)$", \
    "AM_hi.png", 1e-3)

# Q2: sin^3(t) and cos^3(t)
func = lambda x : (sin(x))**3
dft(func, -4*pi, 4*pi, 128, 10, \
    r"$|Y|$", r"$\angle Y$", r"Spectrum of $\sin^{3}(t)$", \
    "sine_cube.png", 1e-3)

func = lambda x : (cos(x))**3
dft(func, -4*pi, 4*pi, 128, 10, \
    r"$|Y|$", r"$\angle Y$", r"Spectrum of $\cos^{3}(t)$", \
    "cos_cube.png", 1e-3)

# Q3: cos(20t + 5cos(t))
func = lambda x : cos(20*x + 5*cos(x))
dft(func, -4*pi, 4*pi, 512, 40, \
    r"$|Y|$", r"$\angle Y$", r"Spectrum of $\cos(20t + 5\cos(t))$", \
    "FM_cos.png", 1e-3)
#'''

# Q4: Gaussian
def gauss(x):
    # Gaussian of mean 0, stddev 1
    return exp((-x**2)/2.0)

def gauss_freq(w):
    # Analog FT of Gaussian outputted by gauss()
    return sqrt(2*pi) * exp(-w**2/2)


T = 0.5*pi
N = 64
err = 999999999
iters = 0
tol = 1e-6

# Stores best [T, w, Y, err]
best = [0,0,0,999999999]

# Iterate over time window sizes, run at least 10 iterations
while err>tol or iters < 10:  
    T*=2
    N*=2

    # time resolution remains same
    t = linspace(-T/2,T/2,N+1)[:-1]

    # frequency resolution increases
    w = linspace(-128,128,N+1)[:-1]

    y = gauss(t)

    # Shift negative time half to front and then do fft, then shift back
    Y = fftshift(fft(ifftshift(y)))*T/N

    err = sum(abs(Y-gauss_freq(w)))/N
    print("Error at iteration " +str(iters)+ ":",err)
    if err < best[3]:
        best = [T,w,Y,err]

    iters += 1
    
print("Error: ",best[3])
print("Time window size (times pi):", best[0]/pi)
Y_0 = gauss_freq(w)

figure()
        
subplot(2,1,1)
plot(w,abs(Y),'r', lw=20)
plot(w,abs(Y_0),'g', lw=5)
xlim([-5,5])
ylabel(r"|Y|",size=16)
legend(['calculated', 'actual'])
title("Best DFT estimate of Gaussian")
grid(True)

subplot(2,1,2)
plot(w,angle(Y),'ro',ms=20, lw=20)
plot(w,angle(Y_0),'go',ms=5, lw=5)
xlim([-5,5])
ylabel(r"$\angle Y$",size=16)
legend(['calculated', 'actual'])
xlabel(r"$k$",size=16)
grid(True)

savefig("best_est.png")
show()