def s24rts(s24,depth,plot,save,sphharm,matrix,lmx=20):
    # THIS PROGRAM READS IN S20RTS AND PLOTS AS MAPS AT DIFFERENT
    # DEPTHS
    import numpy as np
    import latychev_out as lt
    import matplotlib.pyplot as plt
    # from mpl_toolkits.basemap import Basemap
    import sharmonics as sh
    import scipy.special as sp
    import numpy.polynomial.legendre as leg
    from scipy.interpolate import interp1d

    # options:
    # s24: string, either 'S20' or 'S40'
    # depth: float, depth slice to be taken (km)
    # plot: integer, either 0 or 1, if 1, figure will be made
    # save: integer, if 1, a png file is saved 
    # sphharm: integer, either 0 or 1, if 1 spherical harmonic
    #                   decomposition will be printed,
    #                   as a 2D matrix (L x M) in the formalism of JMX.
    # matrix: integer, either 0 or number of lats to be sampled.
    #                  if nonzero, 2D map will be formed from -pi to pi lon.
    # lmx: integer, maximum spherical harmonic degree (must be <=20 for 
    #               S20 and <=40 for S40.
    #               optional. default is lmx=20

    # 1: read in model:
    dir='./'
    f=open(dir+s24+"RTS.dat",'r')
#    f=open(s24+"RTS.dat",'r')
    data=f.read()
    data_split=data.split()
    data_num=[float(i) for i in data_split]
    f.close()
    #if (s24=='S20'):
    #    lmx=20
    #elif (s24=='S40'):
    #    lmx=40
    #else:
    #    print 'wrong file name!'
    #    return None
    ndmx=21
    dv_sph=np.zeros((ndmx,lmx+1,2*lmx+1))
    counter=0
    for k in range(ndmx):
        for l in range(lmx+1):
            for m in range(2*l+1):
                dv_sph[k,l,m]=data_num[counter]
                counter=counter+1

    radchoice=6371.-depth


    # 2: interpolate radially, using cubic spline
    # radial nodes:
    radius=np.array([3480.0000, 3786.2178, 4064.5637, 4317.5884, 4547.5996, 4756.6738,\
                         4946.7188, 5119.4810, 5276.5088, 5419.2646, 5549.0225, 5666.9727,\
                         5774.1899, 5871.6626, 5960.2510, 6040.7852, 6113.9971, 6180.5459,\
                         6241.0327, 6296.0171, 6346.0000])
    dv_slice=np.zeros((lmx+1,2*lmx+1))
    for l in range(lmx+1):
        for m in range(2*l+1):
            f=interp1d(radius,dv_sph[::-1,l,m],kind='cubic')
            dv_slice[l,m]=f(radchoice)

    # 2: set up dimensions of 2D map of model
    # for some reason I'm 90 degrees off in longitude!

    if (matrix==0):
        Nlat=5*lmx # just resolved enough for each model
    else:
        Nlat=matrix

    Nlon=2*Nlat
    [colats,dx]=leg.leggauss(Nlat)
    colats=np.arccos(colats)
    lons=np.linspace(-np.pi,np.pi,Nlon+1) 
    lons=lons[0:-1]
    dv=np.zeros((Nlat,Nlon))


    # 3: calculate ylms at grids
    for l in range(lmx+1):
        #m=0
        ylm=sh.real_sph(colats,lons,l,0)
        dv=dv+dv_slice[l,0]*ylm
        counter=0
        for m in range(1,l+1):
            counter=counter+1
            ylm=sh.real_sph(colats,lons,l,-m)
            dv=dv+dv_slice[l,counter]*ylm
            counter=counter+1
            ylm=sh.real_sph(colats,lons,l,m)
            dv=dv+dv_slice[l,counter]*ylm

    # 4: OUTPUTS
    if (plot==1):
        plt.figure()
        x=np.linspace(0,360,Nlon)
        y=np.linspace(-90.0,90.0,Nlat)
        [XX,YY]=np.meshgrid(x,y)
        m = Basemap(projection='robin',lon_0=180.0,lat_0=0.0,llcrnrlon=0.0,\
                        llcrnrlat=-90,urcrnrlon=360.,urcrnrlat=90.0)
        xx,yy=m(XX,YY)

        m.pcolor(xx,yy,-dv,vmin=-.02,vmax=.02)
        m.drawcoastlines()
        cbar=plt.colorbar(orientation='horizontal',fraction=0.07)
        plt.title('depth '+str(depth)+'km')
        if (save==1):
            plt.savefig(s24+str(depth)+'.png')
        else:
            plt.show()
    
    if (sphharm==1):
        Ulm=sh.fwsh(dv,colats,lmx,dx)
    else: 
        Ulm=0.

    if (matrix==0):
        dv=0.    


    if ((sphharm==1)or(matrix!=0)):
        return Ulm,dv,colats
