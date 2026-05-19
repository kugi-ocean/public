'''2次元浅水方程式の図示: アニメーション'''
import numpy as np
import matplotlib.pyplot as plt

#- ファイル読み込み
f=open('shallow2d.dat','rb')
(im,jm,nm) = np.fromfile(f,'i8',3)
(dx_m,dy_m,dt_sec) = np.fromfile(f,'d',3)
dat = np.fromfile(f,'d',(im+1)*(jm+1)*(nm+1)*3).reshape((im+1,jm+1,3,nm+1),order='F')
f.close()

#- 海面高度
eta = dat[:,:,0,:]
x = (np.arange(im)+0.5)*dx_m * 1.e-3
y = (np.arange(jm)+0.5)*dy_m * 1.e-3

#- 流速 (ベクトルの長さは同じにして、色で流速を示す)
dec = 5  # decimate (間引く間隔)
u = dat[::dec,::dec,1,:]
v = dat[::dec,::dec,2,:]
u_abs=np.sqrt( pow(u,2) + pow(v,2) )
u = u / u_abs
v = v / u_abs
xu = np.arange(im+1)[::dec] * dx_m * 1.e-3
yu = np.arange(jm+1)[::dec] * dy_m * 1.e-3

#- 時間
t = np.arange(nm+1) * dt_sec

#- キャンパス
fig, ax = plt.subplots()

for i in range(0,nm+1,10):  # 3つ目の引数が描画の間隔
    ax.cla()
    cs = ax.contour(x,y,eta[:-1,:-1,i].T)
    #- index order of eta is reversed, so eta.T (transpose) are used.
    Q = ax.quiver(xu, yu, u[:,:,i].T, v[:,:,i].T, u_abs[:,:,i].T, cmap='jet', pivot='mid' )
#    plt.colorbar(Q, label='Velocity [m/s]', shrink=0.6, ax=ax) うまくいかない
    ax.clabel(cs)
    ax.set_title( f'SSH [m] t = {t[i]}',fontsize=10)
    ax.set_xlabel('X [km]')
    ax.set_ylabel('Y [km]')
    plt.pause(0.1)



