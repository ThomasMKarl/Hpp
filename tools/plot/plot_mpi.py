from __future__ import unicode_literals
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.use("pgf")
pgf_with_rc_fonts = {
    "font.family": "serif",
    "font.serif": [],                   # use latex default serif font
    "font.sans-serif": ["DejaVu Sans"], # use a specific sans-serif font
}
mpl.rcParams.update(pgf_with_rc_fonts)
from matplotlib2tikz import save as tikz_save 


font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}

mpl.rc('font', **font)
plt.figure(figsize=(6.5,6))

data = np.loadtxt('ising1.dat')
y1 = data[0:6,1]
sig_y1 = data[0:6,2]
x = data[0:6,0]
fit1 = 6.53*10**-5*x**1.96614

data = np.loadtxt('ising4.dat')
y4 = data[0:6,1]
sig_y4 = data[0:6,2]
fit4 = 3.03*10**-5*x**1.85421

data = np.loadtxt('ising9.dat')
y9 = data[0:6,1]
sig_y9 = data[0:6,2]
fit9 = 4.96*10**-5*x**1.66343


plt.yscale('log', nonposy='clip')
plt.xscale('log', nonposy='clip')
plt.axis([26, 210, 0.006, 2.5])

plt.plot(x, y1, 'o', markerfacecolor='white',markeredgecolor='red', label='1 core')
plt.plot(x, y1, 'o', markerfacecolor='red',markeredgecolor='none',alpha=0.5)
plt.plot(x, y4, 'o', markerfacecolor='white',markeredgecolor='blue', label='4 cores')
plt.plot(x, y4, 'o', markerfacecolor='blue',markeredgecolor='none',alpha=0.5)
plt.plot(x, y9, 'o', markerfacecolor='white',markeredgecolor='green', label='6 cores')
plt.plot(x, y9, 'o', markerfacecolor='green',markeredgecolor='none',alpha=0.5)
#plt.errorbar(x, y1, yerr=sig_y1, fmt='none', elinewidth='2', label='1 core')
#plt.errorbar(x, y4, yerr=sig_y4, fmt='none', elinewidth='2', label='1 core')
#plt.errorbar(x, y9, yerr=sig_y9, fmt='none', elinewidth='2', label='1 core')
plt.plot(x, fit1, color='red',   label='fit(x)=a*x**b')
plt.plot(x, fit4, color='blue',  label='fit(x)=a*x**b')
plt.plot(x, fit9, color='green', label='fit(x)=a*x**b')

plt.ylabel('time / seconds')
plt.xlabel('grid size')
plt.legend(loc='lower right', numpoints=1)

#plt.show()

plt.savefig('ising.pgf', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='pgf', transparent=True, bbox_inches=None, pad_inches=0.1, frameon=None)
