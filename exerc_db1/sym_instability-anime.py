'''浮力振動の図示: アニメーション'''
import numpy as np
import matplotlib.pyplot as plt

f=open('sym.dat','rb')
nm = np.fromfile(f,'i8',1)[0]
h = np.fromfile(f,'d',1)[0]
dat = np.fromfile(f,'d',(nm+1)*4).reshape((4,nm+1),order='F')
y = dat[2]
z = dat[3]
t = np.arange(nm+1)*h
f.close()

xmin=-10.
xmax=10.
xc=np.linspace(xmin,xmax,20)
yc, zc = np.meshgrid(xc,xc)
ac = (-yc + zc)*1.e-6
bc = (-1.2*yc + zc)*1.e-6

fig, ax = plt.subplots()
plt.axis('scaled')

for i in range(nm+1):
    ax.cla()
    ax.set_ylim(xmin, xmax)
    ax.set_xlim(xmin, xmax)
    ax.contour(yc,zc,ac,colors='orange')
    ax.contour(yc,zc,bc,colors='yellow')
    ax.plot(y[i],z[i],'.')
    ax.set_title( f't = {t[i]}',loc='right',fontsize=10)
    plt.pause(0.1)



