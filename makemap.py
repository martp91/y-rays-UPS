import numpy as np
import pylab as plt

def makemap(l,b,values,ps=1,x=360,y=180):
    """Make map x*y with pixelsize=ps from 3 arrays (1d) with longitude latitude
    and the values at l,b
    """
    ny = int(y/ps)
    nx = int(x/ps)
    _map = np.zeros((ny,nx))
    l = l+180
    l = (x)*l/(ps*360.)
    b = b + 90
    b = (y)*b/(ps*180)
    l = l.astype(int)
    b = b.astype(int)
    np.add.at(_map,(b,l),values)
    return _map
    
    
    
#for iy in range(len(_map)):
#        mask_b = (b<=(90-iy*ps))&(b>(90-(iy+1)*ps))
#        for ix in range(len(_map[iy,:])):
#            mask_l = (l>=(-180+ix*ps))&(l<(-180+(ix+1)*ps))
#            mask_lb = mask_b&mask_l
#            val_xy = np.sum(values[mask_lb])
#            _map[iy,ix]= val_xy
#
#
#n=1000
#l = (np.random.rand(n)-0.5)*360
#b = (np.random.rand(n)-0.5)*180
#flux = np.random.rand(n)
#
#makemap(l,b,flux,0.5,360,180)
#
#plt.show()
