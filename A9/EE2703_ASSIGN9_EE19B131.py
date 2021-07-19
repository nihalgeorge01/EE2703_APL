'''
Assignment 9: Spectra of Non-Periodic Signals

Author: Nihal John George (EE19B131)
'''

from pylab import *
import seaborn as sns 
import mpl_toolkits.mplot3d.axes3d as p3

def get_dft(f, T_L, T_R, N, title_str, plot_flag = True, limx = 10, wnd_func = None, img_path = None, func_mode=True):
    
    t = linspace(T_L, T_R, N+1)[:-1]
    dt = t[1] - t[0]
    fmax = 1/dt
    n = arange(N)
    flt_N = float(N)

    if func_mode:
        y = f(t)
    else:
        y = f
    
    if wnd_func != None:
        wnd = wnd_func(n/flt_N)
        y = y*wnd

    y[0] = 0 # sample corresponding to -tmax should be set to 0
    y = fftshift(y)
    Y = fftshift(fft(y))/flt_N
    w = linspace(-pi*fmax, pi*fmax,  N+1)[:-1]

    if plot_flag:
        figure()
        subplot(2, 1, 1)
        plot(w, abs(Y), lw = 2)
        xlim([-limx, limx])
        ylabel(r"$|Y|$", size = 16)
        title(title_str)
        grid(True)

        subplot(2, 1, 2)
        plot(w, angle(Y), 'ro', lw = 2)
        xlim([-limx, limx])
        ylabel(r"$\angle Y$", size = 16)
        xlabel(r"$\omega$", size = 16)
        grid(True)

        savefig(img_path)
        show()

    return Y

FIG_CT = 1
#'''
# Q1: Examples

t = linspace(-pi, pi, 65)[:-1]
dt = t[1]-t[0]
fmax = 1/dt
y = sin(sqrt(2)*t)
y[0] = 0 # the sample corresponding to -tmax should be set zero
y = fftshift(y) # make y start with y(t = 0)
Y = fftshift(fft(y))/64.0
w = linspace(-pi*fmax, pi*fmax, 65)[:-1]

figure(FIG_CT)

subplot(2, 1, 1)
plot(w, abs(Y), lw = 2)
xlim([-10, 10])
ylabel(r"$|Y|$", size = 16)
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)$")
grid(True)

subplot(2, 1, 2)
plot(w, angle(Y), 'ro', lw = 2)
xlim([-10, 10])
ylabel(r"Phase of $Y$", size = 16)
xlabel(r"$\omega$", size = 16)
grid(True)

savefig("fig10-1.png")
show()
FIG_CT += 1

t1 = linspace(-pi, pi, 65)[:-1]
t2 = linspace(-3*pi, -pi, 65)[:-1]
t3 = linspace(pi, 3*pi, 65)[:-1]
# y = sin(sqrt(2)*t)

figure(FIG_CT)

plot(t1, sin(sqrt(2)*t1), 'b', lw = 2)
plot(t2, sin(sqrt(2)*t2), 'r', lw = 2)
plot(t3, sin(sqrt(2)*t3), 'r', lw = 2)
ylabel(r"$y$", size = 16)
xlabel(r"$t$", size = 16)
title(r"$\sin\left(\sqrt{2}t\right)$")
grid(True)

savefig("fig10-2.png")
show()
FIG_CT += 1


t1 = linspace(-pi, pi, 65)[:-1]
t2 = linspace(-3*pi, -pi, 65)[:-1]
t3 = linspace(pi, 3*pi, 65)[:-1]
y = sin(sqrt(2)*t1)

figure(FIG_CT)

plot(t1, y, 'bo', lw = 2)
plot(t2, y, 'ro', lw = 2)
plot(t3, y, 'ro', lw = 2)
ylabel(r"$y$", size = 16)
xlabel(r"$t$", size = 16)
title(r"$\sin\left(\sqrt{2}t\right)$ with $t$ wrapping every $2\pi$ ")
grid(True)

savefig("fig10-3.png")
show()
FIG_CT += 1


t = linspace(-pi, pi, 65)[:-1]
dt = t[1]-t[0]
fmax = 1/dt
y = t
y[0] = 0 # the sample corresponding to -tmax should be set zero
y = fftshift(y) # make y start with y(t = 0)
Y = fftshift(fft(y))/64.0
w = linspace(-pi*fmax, pi*fmax, 65)[:-1]

figure(FIG_CT)

semilogx(abs(w), 20*log10(abs(Y)), lw = 2)
xlim([1, 10])
ylim([-20, 0])
xticks([1, 2, 5, 10], ["1", "2", "5", "10"], size = 16)
ylabel(r"$|Y|$ (dB)", size = 16)
title(r"Spectrum of a digital ramp")
xlabel(r"$\omega$", size = 16)
grid(True)

savefig("fig10-4.png")
show()
FIG_CT += 1

t1 = linspace(-pi, pi, 65)[:-1]
t2 = linspace(-3*pi, -pi, 65)[:-1]
t3 = linspace(pi, 3*pi, 65)[:-1]
n = arange(64)
wnd = fftshift(0.54+0.46*cos(2*pi*n/63))
y = sin(sqrt(2)*t1)*wnd

figure(FIG_CT)

plot(t1, y, 'bo', lw = 2)
plot(t2, y, 'ro', lw = 2)
plot(t3, y, 'ro', lw = 2)
ylabel(r"$y$", size = 16)
xlabel(r"$t$", size = 16)
title(r"$\sin\left(\sqrt{2}t\right)\times w(t)$ with $t$ wrapping every $2\pi$ ")
grid(True)

savefig("fig10-5.png")
show()
FIG_CT += 1

f = lambda x: sin(sqrt(2)*x)
wnd_func = lambda x: fftshift(0.54+0.46*cos(2*pi*x))
Y = get_dft(f, -pi, pi, 64, \
    r"Spectrum of $\sin\left(\sqrt{2}t\right)\times w(t)$", \
    True, 8, wnd_func, "fig10-6.png")

f = lambda x: sin(sqrt(2)*x)
wnd_func = lambda x: fftshift(0.54+0.46*cos(2*pi*x))
Y = get_dft(f, -4*pi, 4*pi, 256, \
    r"Spectrum of $\sin\left(\sqrt{2}t\right)\times w(t)$", \
    True, 8, wnd_func, "fig10-7.png")

# Q2-a: cos^3(wt) without Hamming
f = lambda x: (cos(0.86*x))**3
Y = get_dft(f, -4*pi, 4*pi, 256, \
    r"Spectrum of $\cos\left(0.86t\right)$", \
    True, 8, None, "Q2_no_ham.png")

# Q2-b: cos^3(wt) with Hamming
f = lambda x: (cos(0.86*x))**3
wnd_func = lambda x: fftshift(0.54+0.46*cos(2*pi*x))
Y = get_dft(f, -4*pi, 4*pi, 256, \
    r"Spectrum of $\cos\left(0.86t\right)$", \
    True, 8, wnd_func, "Q2_ham.png")
#'''

# Q3, Q4: cos(wt + d) spectrum estimation, with and without AWGN
#seed(1)
N = 128
wnd_func = lambda x: fftshift(0.54+0.46*cos(2*pi*x))
t = linspace(-pi, pi, N+1)[:-1]
dt = t[1] - t[0]
fmax = 1/dt

w0 = rand(1) + 0.5
d = (rand(1)-0.5)*2*pi

y_clean = cos(w0 * t + d)
y_noisy = y_clean + 0.1*randn(N)
w = linspace(-pi*fmax, pi*fmax, N+1)[:-1]

print("Actual w0, d =", w0, d)

Y_clean = get_dft(y_clean, -pi, pi, N, \
    "Spectrum of cos(%.2f t + %.2f)" % (w0, d), \
    True, 10, wnd_func, "Q3.png", False)
Y_noisy = get_dft(y_noisy, -pi, pi, N, \
    "Spectrum of cos(%.2f t + %.2f)" % (w0, d), \
    True, 10, wnd_func, "Q4.png", False)

# Set low magnitudes to 0 for better estimation
abs_Y_clean = abs(Y_clean)
abs_Y_clean[abs_Y_clean<1e-3] = 0

# Set low magnitudes to 0 for better estimation
abs_Y_noisy = abs(Y_noisy)
abs_Y_noisy[abs_Y_noisy<5e-2] = 0

# Get average frequency weighted by magnitude
w0_est_clean = average(w[N//2:N//2+5], weights=abs_Y_clean[N//2:N//2+5])
w0_est_noisy = average(w[N//2:N//2+5], weights=abs_Y_noisy[N//2:N//2+5])

# Get average phase weighted by magnitude
d_est_clean = (average(angle(Y_clean[N//2+1:]), weights=abs_Y_clean[N//2+1:]) \
                - average(angle(Y_clean[:N//2]), weights=abs_Y_clean[:N//2]))/2
d_est_noisy = (average(angle(Y_noisy[N//2+1:]), weights=abs_Y_noisy[N//2+1:]) \
                - average(angle(Y_noisy[:N//2]), weights=abs_Y_noisy[:N//2]))/2

print("Estimated w0 (clean, noisy):", w0_est_clean, w0_est_noisy)
print("Estimated d  (clean, noisy):", d_est_clean, d_est_noisy)


#'''
# Q5: Chirp
N = 1024
f = lambda x: cos(16*x*(1.5 + x/(2*pi)))
Y = get_dft(f, -pi, pi, N, \
    r"Spectrum of $cos(16(1.5+\frac{t}{2\pi})t)$", \
    True, 64, wnd_func, "Q5_wind.png")

# Q6: Time-Frequency Plot
t = linspace(-pi, pi, N+1)[:-1]
dt = t[1] - t[0]
fmax = 1/dt
blocks = split(t, 16)       # split into 16 blocks, 1024/64 = 16
edges = [[blk[0], blk[63]] for blk in blocks]   # Get block edges

# Find DFT of each block
Y_arr = array([get_dft(f, edge[0], edge[1], 64, \
    r"Spectrum of $cos(16(1.5+\frac{t}{2\pi})t)$", \
    False, 64, wnd_func, "") for edge in edges])

t = t[::64]
w=np.linspace(-fmax*np.pi,fmax*np.pi,65)[:-1]

t_g,w_g = meshgrid(t,w)

fig = figure()
ax = p3.Axes3D(fig)
surf = ax.plot_surface(w_g,t_g, abs(Y_arr.T), cmap=cm.jet)      # Reverse the columns of phi.T to get correct coordinates on the plot
# Formatting to truncate decimal numbers
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.2f'))
xlabel(r"Frequency $\rightarrow$")
ylabel(r"Time $\rightarrow$")
ax.set_zlabel(r"Magnitude $\rightarrow$")
ax.set_title("Surface Plot of Magnitude")
fig.colorbar(surf, shrink=0.5)
savefig('surface_wind.png')
show()

# Spectrogram
figure()
ax = sns.heatmap(abs(Y_arr.T), xticklabels=round_(t,2), yticklabels=round_(w,2), cmap='jet', fmt=".2%")
ax.set_title('Heat Map')
xlabel(r"time $\rightarrow$")
ylabel(r"frequency $\rightarrow$")
savefig('temp_wind.png')
show()
#'''
