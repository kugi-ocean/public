'''浮力振動の図示: t-zプロット'''
import numpy as np
import matplotlib.pyplot as plt

f=open('buoyancy.dat','rb')
nm = np.fromfile(f,'i8',1)[0]
h = np.fromfile(f,'d',1)[0]
z = np.fromfile(f,'d',(nm+1)*2).reshape((2,nm+1),order='F')[0]
t = np.arange(nm+1)*h
f.close() 

fig, ax = plt.subplots()
ax.plot(t,z)
ax.set_xlabel('t')
ax.set_ylabel('z')
plt.show()





