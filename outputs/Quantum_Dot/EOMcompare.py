import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')

filename1 = 'SarahQDdata.dat'
filename2 = 'PA_final.dat'

cases1 = {('6','1.0'):0,('6','0.28'):1,('6','0.1'):2,('12','1.0'):3,('12','0.28'):4,('12','0.1'):5,('20','1.0'):6,('20','0.28'):7,('20','0.1'):8,('30','1.0'):9,('30','0.28'):10,('30','0.1'):11,('42','1.0'):12,('42','0.28'):13,('42','0.1'):14}
cases2 = {('2','1.00'):0,('2','0.28'):1,('2','0.10'):2,('3','1.00'):3,('3','0.28'):4,('3','0.10'):5,('4','1.00'):6,('4','0.28'):7,('4','0.10'):8,('5','1.00'):9,('5','0.28'):10,('5','0.10'):11,('6','1.00'):12,('6','0.28'):13,('6','0.10'):14}

shells = {'2':'6','3':'12','4':'20','5':'30','6':'42'}


x1 = [[] for i in range(15)]
y1 = [[] for i in range(15)]
x2 = [[] for i in range(15)]
y2 = [[] for i in range(15)]

ind = -1
for line in open(filename1):
    if(line == '\n'):
        continue
    elif(line.find('####') > -1):
        line = line.split()
        ind = cases1.get((line[1],line[2]), -1)
        #print(line,'ind = ',ind)
    else:
        if(ind == -1):
            continue
        line = line.split()
        #print(line,len(line))
        if(len(line) < 4):
            continue
        else:
            #print(line[0],line[3])
            x1[ind].append(float(line[0]))
            y1[ind].append(float(line[3]))


with open(filename2) as f2:
    data2 = f2.read()
data2 = data2.split('\n')
ind = -1
for num in range(2,len(data2)):
    if((num - 2)%3 == 0):
        line = data2[num].split()
        ind = cases2.get((line[1],line[4]), -1)
        print(line)
        print(line[0],line[4],'ind = ',ind)
        print(line[1],line[5])
        if(ind == -1):
            continue
        x2[ind].append(float(line[0]))
        y2[ind].append(float(line[5]))


#print(len(x1[0]),len(y1[0]),len(x2[0]),len(y2[0]))
#print(' '.join(map(str,x1[0])))
#print(' '.join(map(str,y1[0])))
#print(' '.join(map(str,x2[0])))
#print(' '.join(map(str,y2[0])))

for num in range(15):
    plt.plot(x1[num],y1[num],x2[num],y2[num])
    plt.savefig('fig'+str(num+1)+'.pdf', format='pdf')
    plt.clf()

#plt.show()
