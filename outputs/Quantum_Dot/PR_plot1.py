import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

majorLocatorX = MultipleLocator(2)
minorLocatorX = MultipleLocator(1)
majorLocatorY = MultipleLocator(0.5)
minorLocatorY = MultipleLocator(0.25)

filename = 'PR_final.dat'
filename2 = 'EOM_IMSRG_FEI_HAM_particle_removed.dat'

x = [[] for i in range(3)]
y = [[] for i in range(3)]

x2 = [[] for i in range(3)]
y2 = [[] for i in range(3)]

shells = int(sys.argv[1])
hw = float(sys.argv[2])

with open(filename) as f:
    data = f.read()
data = data.split('\n')

for num in range(2,len(data)):
    line = data[num].split()
    if(int(line[1]) == shells and float(line[4]) == hw):
        ind = int(line[2])
        x[ind].append(float(line[0]))
        y[ind].append(float(line[6]))


with open(filename2) as f2:
    data2 = f2.read()
data2 = data2.split('\n')
    
for num in range(2,len(data2)):
    line = data2[num].split()
    if(int(line[1]) == shells and float(line[4]) == hw):
        ind = int(line[2])
        x2[ind].append(float(line[0]))
        y2[ind].append(float(line[6]))
        
        
plt.rc('font', family='serif')

fig, ax = plt.subplots()
ax.plot(x[0], y[0], '-', marker='o', color='k', label='$\mathrm{CCSD: m_{l}=0}$')
ax.plot(x[1], y[1], '--', marker='s', color='k', label='$\mathrm{CCSD: m_{l}=1}$')
ax.plot(x[2], y[2], ':', marker='^', color='k', label='$\mathrm{CCSD: m_{l}=2}$')
ax.plot(x2[0], y2[0], '-', marker='o', color='r', label='$\mathrm{IM-SRG: m_{l}=0}$')
ax.plot(x2[1], y2[1], '--', marker='s', color='r', label='$\mathrm{IM-SRG: m_{l}=1}$')
ax.plot(x2[2], y2[2], ':', marker='^', color='r', label='$\mathrm{IM-SRG: m_{l}=2}$')
ax.legend(loc='lower right')

plt.xlabel(r'$\mathrm{Total Shells}$', fontsize=15)
plt.ylabel(r'$\mathrm{\Delta E_{k}}$', fontsize=15)
plt.axis([5, 21, -13.75, -12])

annotation_string = r'$\mathrm{\hbar\omega=1.00}$'
annotation_string += '\n'
annotation_string += r'$\mathrm{20\ particles}$'
ax.annotate(annotation_string, fontsize=15, xy=(6, -12.25))

ax.xaxis.set_major_locator(majorLocatorX)
ax.xaxis.set_minor_locator(minorLocatorX)
ax.yaxis.set_major_locator(majorLocatorY)
ax.yaxis.set_minor_locator(minorLocatorY)

plt.savefig('PRfig_'+str(shells)+'_'+str(hw)+'.pdf', format='pdf', bbox_inches='tight')
plt.show()
