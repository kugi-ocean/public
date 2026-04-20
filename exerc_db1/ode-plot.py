'''常微分方程式の計算解の図示: t-xプロット'''
import numpy as np
import matplotlib.pyplot as plt

f=open('ode.dat','rb')
nm = np.fromfile(f,'i8',1)[0]
h = np.fromfile(f,'d',1)[0]
x = np.fromfile(f,'d',nm*2).reshape((2,nm),order='F')
t = np.arange(nm)*h
f.close() 

fig, ax = plt.subplots()
ax.plot(t,x[1],label='true')
ax.plot(t,x[0],label='calc')
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.legend()
plt.show()





