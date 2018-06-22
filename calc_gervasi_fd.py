import numpy as np

const_c = 3e8
const_k = 1.38e-23

gervasi = open('gervasi.txt','r')

freq=[];temp=[]
while True:
	line = gervasi.readline()
	if not line: break

	if line[0]!='#':
		items = line.split()

		freq.append(float(items[0])) # MHz
		temp.append(float(items[1])) # mK 

freq = np.array(freq)
temp = np.array(temp)

lam = const_c/(freq*1e6)

fd = 1e26*2*const_k*temp*1e-3/lam**2

for i in range(0,len(fd)):
	print freq[i],fd[i]