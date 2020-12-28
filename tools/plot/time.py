from __future__ import unicode_literals
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.use("pgf")
plt.figure(figsize=(6,6))

data1 = np.loadtxt('time1.dat')
data4 = np.loadtxt('time4.dat')
y1 = data1[:,1]
x  = data1[:,0]*data1[:,0]*data1[:,0]
y2 = data4[:,1]

#plt.yscale('log', nonposy='clip')
#plt.xscale('log', nonposy='clip')
#plt.axis([26, 210, 0.006, 2.5])

plt.plot(x, y1, 'o', markerfacecolor='white',markeredgecolor='red', label='1 core')
plt.plot(x, y1, 'o', markerfacecolor='red',markeredgecolor='none',alpha=0.5)
plt.plot(x, y2, 'o', markerfacecolor='white',markeredgecolor='blue', label='6 cores')
plt.plot(x, y2, 'o', markerfacecolor='blue',markeredgecolor='none',alpha=0.5)

#plt.plot(x, fit, color='red',   label='fit(x)=a*x**b')
#plt.plot(x, fit2, color='blue',  label='fit(x)=a*x**b')

plt.ylabel('Time / sec')
plt.xlabel('Volume')
plt.legend(loc='upper left', numpoints=1)

plt.savefig('time.pgf', dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='pgf', transparent=True, bbox_inches=None, pad_inches=0.1, frameon=None)

plt.show()
