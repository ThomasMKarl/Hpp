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


font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}

mpl.rc('font', **font)
plt.figure(figsize=(6,6))

data = np.loadtxt('ising.dat')
y = data[:,1]
sig_y = data[:,2]
x = data[:,0]

#plt.yscale('log', nonposy='clip')
#plt.xscale('log', nonposy='clip')
#plt.axis([26, 210, 0.006, 2.5])

plt.plot(x, y, 'o', markerfacecolor='white',markeredgecolor='red', label='')
plt.plot(x, y, 'o', markerfacecolor='red',markeredgecolor='none',alpha=0.5)
#plt.plot(x, y2, 'o', markerfacecolor='white',markeredgecolor='blue', label='6 cores')
#plt.plot(x, y2, 'o', markerfacecolor='blue',markeredgecolor='none',alpha=0.5)

#plt.plot(x, fit, color='red',   label='fit(x)=a*x**b')
#plt.plot(x, fit2, color='blue',  label='fit(x)=a*x**b')

plt.ylabel('Energy per Spin / J')
plt.xlabel('Temperature / $J/k_b$')
plt.legend(loc='upper left', numpoints=1)

plt.savefig('ising.pgf', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='pgf', transparent=True, bbox_inches=None, pad_inches=0.1, frameon=None)

plt.show()
