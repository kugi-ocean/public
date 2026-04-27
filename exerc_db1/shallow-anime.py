'''1次元浅水方程式の図示: アニメーション'''
import numpy as np
import matplotlib.pyplot as plt

f=open('shallow.dat','rb')
(mm,nm) = np.fromfile(f,'i8',2)
(dx_m,h) = np.fromfile(f,'d',2)
dat = np.fromfile(f,'d',(mm+1)*(nm+1)*3).reshape((mm+1,3,nm+1),order='F')
#dat = np.fromfile(f,'d',(mm+1)*(nm+1)*2).reshape((mm+1,2,nm+1),order='F')
eta = dat[:,0,:]
x = (np.arange(mm)+0.5)*dx_m
t = np.arange(nm+1)*h
f.close()

fig, ax = plt.subplots()

for i in range(nm+1):
    ax.cla()
#    ax.set_ylim(-0.01, 0.01)
    ax.plot(x,eta[:-1,i])
    ax.set_title( f't = {t[i]}',loc='right',fontsize=10)
    plt.pause(0.1)



