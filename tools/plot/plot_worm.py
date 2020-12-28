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
        'size'   : 16}

mpl.rc('font', **font)
plt.figure(figsize=(6,6.2))

data  = np.loadtxt('output_worm.dat')
temp  = data[0:100,0]
en    = data[0:100,1]
sig_e = data[0:100,2]
m     = data[0:100,3]
sig_m = data[0:100,4]



plt.plot(temp, en, 'o', markerfacecolor='white',markeredgecolor='red')
plt.plot(temp, en, 'o', markerfacecolor='red',markeredgecolor='none',alpha=0.5)
#plt.axis([0, 210, 0, 1.5,])
#plt.title(r'$\mu=0.5,\ \beta=1$')
plt.ylabel('Energy per Volume / J')
plt.xlabel('Temperature / J/kb')

#plt.show()

plt.savefig('worme.pgf', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='pgf', transparent=True, bbox_inches=None, pad_inches=0.1, frameon=None)


#plt.errorbar(temp, m, yerr=sig_m, fmt='none', elinewidth='2')
#plt.yscale('log', nonposy='clip')

#plt.ylabel('Magnetic Susceptibility')
#plt.xlabel('Temperature / ')

#plt.show()

#plt.savefig('wormm.pgf', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='pgf', transparent=True, bbox_inches=None, pad_inches=0.1, frameon=None)
