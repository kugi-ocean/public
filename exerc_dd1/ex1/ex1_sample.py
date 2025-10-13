import numpy as np
import matplotlib.pyplot as plt

devtyp='X'
#devtyp='PDF'

# Set parameters

rdir='./'
expid = 'ex1-case1'

nbgn =
nend =
nskp = 1


## Set coordinates

km =
dz =
z = (np.arange(km+2)-0.5 - km) * dz 
dtyp = np.dtype( [('time','<d'), \
                  ('u','<'+str(km+2)+'d'), \
                  ('v','<'+str(km+2)+'d')] )

## Main loop

for n in range(nbgn,nend+1,nskp):

    ### read data
    
    fname = rdir + expid + f'.n{n:06}'
    fp=open(fname,'rb')
    chunk = np.fromfile(fp, dtype=dtyp)
    fp.close()

    time = chunk['time'][0]
    u = chunk['u'][0]
    v = chunk['v'][0]

    ### Make figure
    fig, (figa, figb) = plt.subplots(figsize=(12,6),ncols=2)

    ### Fig1: Z Profile
    figa...

    ### Fig2: Hodograph
    
    figb...


    fig.tight_layout()

    if devtyp == 'X':
        fig.show()
    elif devtyp == 'PDF':
        fig.savefig(edir+expid+f'.n{n:03}'+'.vel.pdf')
        plt.close()
