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

xmax = np.abs([x,y]).max()  #- 図示する範囲を最大値に合わせる
xmin = -xmax

fig, ax = plt.subplots()
plt.axis('scaled')
ax.set_ylim(xmin, xmax)
ax.set_xlim(xmin, xmax)
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')

for i in range(nm+1):
    ax.plot(x[i],y[i],'.',color='black')
    ax.set_title( f't = {t[i]} [sec]',loc='right',fontsize=10)
    plt.pause(0.1)

plt.show()

