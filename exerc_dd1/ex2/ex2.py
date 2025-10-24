'''海面高度と流速ベクトルを描く'''

import numpy as np
import matplotlib.pyplot as plt
import time as tm

devtyp='png'  # x/png/PDF

# set parameters

rdir='./'
expid = 'ex2-case1'

nbgn = 36
nend = 36
nskp = 1

im = 150
km = 15
dx_km = 100.e0  # [km]
dz_m = 100.e0   # [m]

## set coordinates

xc = (np.arange(im+2)-0.5) * dx_km
zc = (np.arange(km+2)-0.5 - km) * dz_m

#- for drawing vectors
xc2 = (np.arange(im)) * dx_km
zc2 = (np.arange(km) - km) * dz_m
ilist=list(range(0,im-1,4))+[im-1]  # ベクトルを間引いて表示する
x2,z2 = np.meshgrid(xc2[ilist], zc2, indexing = 'ij') 

dtyp = np.dtype( [('time','<d'),\
                  ('u','<'+str((im+2)*(km+2))+'d'),\
                  ('v','<'+str((im+2)*(km+2))+'d'),\
                  ('e','<'+str((im+2))+'d'),\
                  ('w','<'+str((im+2)*(km+2))+'d'),\
                  ('p','<'+str((im+2)*(km+2))+'d')] )

## Main Loop

for n in range(nbgn,nend+1,nskp):

    ### read data
    fname = rdir + expid + f'.n{n:06}'
    fp=open(fname,'rb')
    chunk = np.fromfile(fp, dtype=dtyp)
    fp.close()

    time = chunk['time'][0]
    u = chunk['u'][0].reshape((im+2,km+2),order="F")
    v = chunk['v'][0].reshape((im+2,km+2),order="F")
    e = chunk['e'][0].reshape((im+2),order="F")
    w = chunk['w'][0].reshape((im+2,km+2),order="F")
    p = chunk['p'][0].reshape((im+2,km+2),order="F")

    ### Make figure
    fig, figa = plt.subplots(figsize=(12,6),ncols=1)

    ### SSH (eta) Profile
    figa.plot(xc,e[:]*10000.)
    figa.grid()

    ### Velocity: ベクトルの長さは同じにして、色で流速を示す
    u2=0.5*(u[0:-2,1:-1]+u[1:-1,1:-1])
    w2=0.5*(w[1:-1,0:-2]+w[1:-1,1:-1])
    u_abs=np.sqrt( pow(u2,2) + pow(w2,2) )
    u2 = u2 / u_abs
    w2 = w2 / u_abs

    Q = figa.quiver(x2,z2,u2[ilist,:],w2[ilist,:],u_abs[ilist,:],cmap='jet',pivot='mid')
    plt.colorbar(Q, label='Velocity [m/s]')

    figa.set_title('Sea Surface Height')
    figa.set_xlabel('X [km]')
    figa.set_ylabel('Z [m]')

    ### Output
    fig.tight_layout()
    if devtyp == 'X':
        fig.show()
    elif devtyp == 'png':
        fig.savefig('temp.png')
        plt.close()
    elif devtyp == 'PDF':
        fig.savefig(edir+expid+f'.n{n:03}'+'.vel.pdf')
        plt.close()
