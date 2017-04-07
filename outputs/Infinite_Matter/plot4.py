#!/usr/bin/python

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

f = open('output_particlesHF.txt', "rb")
data = f.read()
data = data.rstrip()
data = data.split('\n')

x0 = [row.split('\t')[2] for row in data]
y0 = [row.split('\t')[4] for row in data]

x1 = [x0[i] for i in range(0, 17)]
y1 = [y0[i] for i in range(0, 17)]
x2 = [x0[i] for i in range(17, 34)]
y2 = [y0[i] for i in range(17, 34)]

for i in range(0,17):
    y1[i] = abs(float(y1[i]))

for i in range(0,17):
    y2[i] = abs(float(y2[i]))

plt.rc('font', family='serif')

fig, ax = plt.subplots()
ax.loglog(x1, y1, '-', marker='o', color='b', label='$\mathrm{T_{kin}}\ (\mathrm{PBC})$')
ax.loglog(x2, y2, '--', marker='s', color='r', label='$\mathrm{T_{kin}}\ (\mathrm{TABC5})$')
ax.legend(loc='lower left')

plt.xlabel(r'$A$', fontsize=15)
plt.ylabel(r'$|1-\mathrm{T_{N}}/\mathrm{T_{inf}}|$', fontsize=15)

ax.annotate(r'$\mathrm{\rho=0.16\ fm^{-3}}$', fontsize=15, xy=(350, 0.04), xytext=(350, 0.04))

plt.savefig('fig4.pdf', format='pdf', bbox_inches='tight')
plt.show()
