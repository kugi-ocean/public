'''慣性振動の図示: x-yプロット'''
import numpy as np
import matplotlib.pyplot as plt

f=open('inertial.dat','rb')
nm = np.fromfile(f,'i8',1)[0]
h = np.fromfile(f,'d',1)[0]
dat = np.fromfile(f,'d',(nm+1)*4).reshape((4,nm+1),order='F')
x = dat[2]
y = dat[3]
t = np.arange(nm+1)*h
f.close() 

fig, ax = plt.subplots()
ax.plot(t,y)
ax.set_xlabel('t[sec]')
ax.set_ylabel('y[m]')
plt.show()





