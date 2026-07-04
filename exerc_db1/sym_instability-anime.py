'''浮力振動の図示: アニメーション'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

f=open('sym.dat','rb')
nm = np.fromfile(f,'i8',1)[0]
h = np.fromfile(f,'d',1)[0]
ay_ps2 = np.fromfile(f,'d',1)[0]
az_ps2 = np.fromfile(f,'d',1)[0]
by_ps2 = np.fromfile(f,'d',1)[0]
bz_ps2 = np.fromfile(f,'d',1)[0]
dat = np.fromfile(f,'d',(nm+1)*4).reshape((4,nm+1),order='F')
y = dat[2]
z = dat[3]
t = np.arange(nm+1)*h
f.close()

xmax = np.abs([y,z]).max()  #- 図示する範囲を最大値に合わせる
xmin = - xmax

xc=np.linspace(xmin,xmax,20)
yc, zc = np.meshgrid(xc,xc)

#- Set ac and bc following exp parameters.
ac = ay_ps2 * yc + az_ps2 * zc
bc = by_ps2 * yc + bz_ps2 * zc

fig, ax = plt.subplots()
plt.axis('scaled')

ax.set_ylim(xmin, xmax)
ax.set_xlim(xmin, xmax)
ax.contour(yc,zc,ac,colors='orange')
ax.contour(yc,zc,bc,colors='yellow')
ax.set_xlabel('y[m]')
ax.set_ylabel('z[m]')

legend_lines = [
    mlines.Line2D([], [], color='orange', label='A(=fM)'),
    mlines.Line2D([], [], color='yellow', label='B(=g log theta)')
]
plt.legend(handles=legend_lines)

for i in range(nm+1):
    ax.plot(y[i],z[i],'.',color='black')
    ax.set_title( f't = {t[i]} [sec]',loc='right',fontsize=10)
    plt.pause(0.1)

plt.show()

