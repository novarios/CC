import sys
import numpy as np
import matplotlib.pyplot as plt
#plt.switch_backend('agg')

filename = 'PR_final.dat'

x = [[] for i in range(3)]
y = [[] for i in range(3)]

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

#for num in range(15):
#    plt.plot(x1[num],y1[num],x2[num],y2[num])
#    plt.savefig('fig'+str(num+1)+'.pdf', format='pdf')
#    plt.clf()

plt.plot(x[0],y[0],x[1],y[1],x[2],y[2])
plt.show()
