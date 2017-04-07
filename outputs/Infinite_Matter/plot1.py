#!/usr/bin/python

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

majorLocatorX = MultipleLocator(0.05)
minorLocatorX = MultipleLocator(0.025)
majorLocatorY = MultipleLocator(2)
minorLocatorY = MultipleLocator(1)

f = open('output_numholes.txt', "rb")
data = f.read()
data = data.rstrip()
data = data.split('\n')

x0 = [row.split('\t')[3] for row in data]
y0 = [row.split('\t')[4] for row in data]

x1 = [x0[i] for i in range(0, 14)]
y1 = [y0[i] for i in range(0, 14)]
x2 = [x0[i] for i in range(14, 28)]
y2 = [y0[i] for i in range(14, 28)]
x3 = [x0[i] for i in range(28, 42)]
y3 = [y0[i] for i in range(28, 42)]
x4 = [x0[i] for i in range(42, 56)]
y4 = [y0[i] for i in range(42, 56)]
x5 = [x0[i] for i in range(56, 70)]
y5 = [y0[i] for i in range(56, 70)]
x6 = [x0[i] for i in range(70, 84)]
y6 = [y0[i] for i in range(70, 84)]

plt.rc('font', family='serif')

fig, ax = plt.subplots()
ax.plot(x1, y1, '-', marker='o', color='k', label='$\mathrm{A=14}$')
ax.plot(x2, y2, '--', marker='s', color='b', label='$\mathrm{A=54}$')
#ax.plot(x3, y3, ':', marker='^', color='r', label='$\mathrm{A=114}$')
ax.plot(x4, y4, '-.', marker='^', color='g', label='$\mathrm{A=186}$')
#ax.plot(x5, y5, ':', marker='>', color='r', label='$\mathrm{A=294}$')
ax.plot(x6, y6, ':', marker='v', color='r', label='$\mathrm{A=358}$')
ax.legend(loc='lower right')

plt.xlabel(r'$\mathrm{\rho\ [fm^{-3}]}$', fontsize=15)
plt.ylabel(r'$\mathrm{E/A\ [MeV]}$', fontsize=15)
plt.axis([0.0, 0.375, 4, 17])

ax.annotate(r'$\mathrm{N_{max}=20}$', fontsize=15, xy=(0.025, 16), xytext=(0.025, 16))

ax.xaxis.set_major_locator(majorLocatorX)
#ax.xaxis.set_major_formatter(majorFormatterX)
ax.xaxis.set_minor_locator(minorLocatorX)
ax.yaxis.set_major_locator(majorLocatorY)
#ax.yaxis.set_major_formatter(majorFormatterY)
ax.yaxis.set_minor_locator(minorLocatorY)

plt.savefig('fig1.pdf', format='pdf', bbox_inches='tight')
plt.show()
