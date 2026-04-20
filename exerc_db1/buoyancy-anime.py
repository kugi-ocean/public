'''浮力振動の図示: アニメーション'''
import numpy as np
import matplotlib.pyplot as plt

f=open('buoyancy.dat','rb')
nm = np.fromfile(f,'i8',1)[0]
h = np.fromfile(f,'d',1)[0]
z = np.fromfile(f,'d',(nm+1)*2).reshape((2,nm+1),order='F')[0]
t = np.arange(nm+1)*h
f.close()

fig, ax = plt.subplots()
plt.axis('scaled')

for i in range(nm+1):
    ax.cla()
    ax.set_ylim(-2., 2.)
    ax.set_xlim(-0.2, 0.2)
    ax.plot(0.0,z[i],'.')
    ax.set_title( f't = {t[i]}',loc='right',fontsize=10)
    plt.pause(0.1)



