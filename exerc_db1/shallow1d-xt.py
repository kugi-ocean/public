'''1次元浅水方程式の図示: x-t プロット'''
import numpy as np
import matplotlib.pyplot as plt

f=open('shallow1d.dat','rb')
(im,nm) = np.fromfile(f,'i8',2)
(dx_m,h) = np.fromfile(f,'d',2)
dat = np.fromfile(f,'d',(im+1)*(nm+1)*2).reshape((im+1,2,nm+1),order='F')
#dat = np.fromfile(f,'d',(im+1)*(nm+1)*3).reshape((im+1,3,nm+1),order='F')  # for (eta,u,v)
eta = dat[:,0,:]
x = (np.arange(im)+0.5)*dx_m
t = np.arange(nm+1)*h
f.close()

ymax = np.abs(eta).max()  #- 図示する範囲を最大値に合わせる

fig, ax = plt.subplots()

cs = ax.contour(x,t,eta[:-1,:].T)
ax.clabel(cs)
ax.set_title( 'eta [m]')
ax.set_xlabel('x [m]')
ax.set_ylabel('t [sec]')

plt.show()
