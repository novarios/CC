filename = 'PA_final.txt'
filename1 = 'PA_final1.txt'
filename2 = 'PR_final1.txt'

with open(filename) as f:
    data = f.read()
data = data.split('\n')

f1 = open(filename1,'w')
f2 = open(filename2,'w')

for num in range(0,len(data)):
    if(num%6 < 3):
        f1.write(data[num])
        f1.write('\n')
    else:
        f2.write(data[num])
        f2.write('\n')

f1.close()
f2.close()
