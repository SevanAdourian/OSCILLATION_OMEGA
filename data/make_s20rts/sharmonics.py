def interpR2GL(X,colats,lons):
    # H Lau Oct 2014
    # interpolates from regular grid to:
    #     - Gauss-Legendre grid in latitude direction
    #     - 2^n nodes in the longitude direction (for fft)
    # use only if your field is smooth enough.
    # minimum degrees
    # inputs:
    # X: 2D array of the field
    # colat: 1D array of colatitude nodes (radians)
    # lons: 1D array of the longitudes (radians)
    # interpkind: type of interpolation

    # outputs:
    # Xnew: interpolated 2D grid
    # newcolats: 1D array of new colatitude nodes (radians)
    # newlons: 1D array of new longitude nodes (radians)

    import numpy as np
    import scipy.special as sp
    import numpy.polynomial.legendre as leg
    from scipy import interpolate
    
    powers=np.arange(4,11) # min=16, max=1024
    idx = (np.abs(2**powers-len(colats))).argmin()
    Nlats=len(colats)
    Nlons=len(lons)
    Nlats_new=2**powers[idx]
    Nlons_new=2*Nlats_new
    [newcolats,dx]=leg.leggauss(Nlats_new)
    newcolats=np.arccos(newcolats)[::-1]
    newlons=np.linspace(0,2.*np.pi,Nlons_new+1)
    newlons=newlons[0:-1]

    x,y=np.meshgrid(lons,colats)
    pts=np.zeros((Nlats*Nlons,2))
    pts[:,0]=np.reshape(x,(Nlats*Nlons))
    pts[:,1]=np.reshape(y,(Nlats*Nlons))
    Xlist=np.reshape(X,(Nlats*Nlons))
    xnew, ynew = np.meshgrid(newlons,newcolats)

    Xnew = interpolate.griddata(pts, Xlist, (xnew, ynew), method='cubic')


    return Xnew,newcolats,newlons,dx

def interpR2GL_v2(X,colats,lons,Nlats_new):
    # H Lau Oct 2014
    # interpolates from regular grid to:
    #     - Gauss-Legendre grid in latitude direction
    #     - 2^n nodes in the longitude direction (for fft)
    # use only if your field is smooth enough.
    # minimum degrees
    # inputs:
    # X: 2D array of the field
    # colat: 1D array of colatitude nodes (radians)
    # lons: 1D array of the longitudes (radians)
    # Nlats: set the number of latitudinal nodes for the interpolation
    #        must be a number of base two

    # outputs:
    # Xnew: interpolated 2D grid
    # newcolats: 1D array of new colatitude nodes (radians)
    # newlons: 1D array of new longitude nodes (radians)

    import numpy as np
    import scipy.special as sp
    import numpy.polynomial.legendre as leg
    from scipy import interpolate
    
    Nlats=len(colats)
    Nlons=len(lons)
    Nlons_new=2*Nlats_new
    [newcolats,dx]=leg.leggauss(Nlats_new)
    newcolats=np.arccos(newcolats)[::-1]
    newlons=np.linspace(0,2.*np.pi,Nlons_new+1)
    newlons=newlons[0:-1]

    x,y=np.meshgrid(lons,colats)
    pts=np.zeros((Nlats*Nlons,2))
    pts[:,0]=np.reshape(x,(Nlats*Nlons))
    pts[:,1]=np.reshape(y,(Nlats*Nlons))
    Xlist=np.reshape(X,(Nlats*Nlons))
    xnew, ynew = np.meshgrid(newlons,newcolats)

    Xnew = interpolate.griddata(pts, Xlist, (xnew, ynew), method='cubic')


    return Xnew,newcolats,newlons,dx




def fwsh(X,colat,lmax,dx):
    # by Harriet Lau; 2 Jul 2013
    # performs the forward spherical harmonic transform of a 2D field.

    # inputs:
    # X: 2D array of the field to transform
    # colat: 1D array of the colatitudes at each point of field X
    #        along vertical dimension (in radians) 
    # lmax: the maximum degree of spherical harmonic to convert
    # dx: the spacing of each point (same size as colat)
    # (this would be in the x=cos(colat) space)

    # outputs:
    # U: 2D array spherical components of lmax+1 by lmax +1
    #    rows are degree l and columns are order m

    # notes:
    # works best when horizontal length of X is a base 2 number
    # and when the vertical length is half that. 
    # e.g. if len(X[0,:])=256, then len(X[0,:])=128
    # also if points colat and dx are Gauss-Legendre roots and weights,
    # respectively

    import numpy as np
    import scipy.special as sp
    import math as ms
    
    # initialize arrays    
    NN=len(X[0,:])
    MM=len(colat)
    G=np.empty((len(colat),lmax+1),dtype=complex)
    G[:,:]=np.NAN
    P=np.empty((MM,lmax+1,lmax+1),dtype=complex)
    P[:,:,:]=np.NAN
    U=np.empty((lmax+1,lmax+1),dtype=complex)
    U[:,:]=np.NAN

    # take FFT accross latitudes
    for i in range(MM):
        FT=np.fft.fft(X[i,:])*(2.*np.pi/NN)
        G[i,0:lmax+1]=FT[0:lmax+1]

    # computer Associated Legendre polynomials
    for l in range(lmax+1):
        for m in range(l+1):
            FACT= np.sqrt((2.*l+1.)*ms.factorial(l-m)/ \
                                        ms.factorial(l+m))
	    P[:,l,m]=sp.lpmv(m,l,np.cos(colat))*FACT
            
            #FACT=(-1.)**m * np.sqrt((2*l+1)*ms.factorial(l-m)/ \
            #                            ms.factorial(l+m))
            #for j in range(MM):
            #    P[j,l,m]=sp.lpmv(m,l,np.cos(colat[j]))*FACT
            
    for l in range(lmax+1):
        P[:,l,:]=P[:,l,:]*G
    
    # integrate in x direction
    for l in range(lmax+1):
        for m in range(l+1):
            U[l,m]= np.sum(P[:,l,m]*dx,axis=0)/(4.*np.pi)
            
    return U


def ivsh(U,colat,lons,lmax):
    # by Harriet Lau; 2 Jul 2013
    # performs the inverse spherical harmonic transform of a spherical
    # harmonic coefficients to produce a global field.

    # inputs:
    # U: 2D array of the spherical harmonic coefficients
    #    of l degrees in rows and m orders in columns
    #    where m = 0,1,...l
    # colat: 1D array of the colatitudes at each point of field X
    #        along vertical dimension (in radians) 
    # lons: longitudes of grid (in radians from 0 to 2pi)
    # lmax: the maximum degree of spherical harmonic to convert
    
    # outputs:
    # X: 2D array of the global field

    import numpy as np
    import scipy.special as sp
    import math as ms

    MM=len(colat)
    NN=len(lons)
    X=np.zeros((MM,NN))
    Y=np.zeros((MM,NN),dtype=complex)
    

    for l in range(lmax+1):
        # first the 0-th order
        for i in range(MM):
            #Y[i,:]=sp.sph_harm(0,l,lons,colat[i])
            m=0
            FACT=np.sqrt((2.*l+1.)*ms.factorial(l-m)/ \
                                        ms.factorial(l+m))            
	    #FACT=np.sqrt(2.)
            Y[i,:]=sp.lpmv(m,l,np.cos(colat[i]))*FACT* \
                np.exp(1.j*m*lons)
        X[:,:]=X[:,:]+np.real(U[l,0])*np.real(Y)

    for l in range(1,lmax+1):
        for m in range(1,l+1):
            for i in range(MM):
            # form spherical harmonics
                FACT= np.sqrt((2.*l+1.)*ms.factorial(l-m)/ \
                                        ms.factorial(l+m))            
                #FACT= np.sqrt(2.)
		Y[i,:]=sp.lpmv(m,l,np.cos(colat[i]))*FACT* \
                    np.exp(1.j*m*lons)

#                Y[i,:]=sp.sph_harm(m,l,lons,colat[i])
            # perform summation    
            X[:,:]=X[:,:]+2.*np.real(U[l,m]*Y)   
                
    return X

def ivsh_rotatelon(U,colat,lons,lmax,dphi):
    # by Harriet Lau; Jul 2016
    # performs the inverse spherical harmonic transform of a spherical
    # harmonic coefficients to produce a global field.
    # but also rotates field by dphi radians.
    # this is useful as some fields output are -180->180 and rotating
    # can be cumbersome.  Easier just to rotate the Ulms.
    # this outputs the 2D field (as ivsh does) but also the
    # updated Ulms that include the rotations.

    # so e.g., just stick in map of -180 to 180 in fwsh.
    # this is not strictly correct as fwsh assumes 0 to 360
    # then use ivsh_rotatelon to produce rotated 2D map (as a check).
    # Should look like 0 to 360.
    # Then we output Ulms adjusted to make such a map.

    # inputs:
    # U: 2D array of the spherical harmonic coefficients
    #    of l degrees in rows and m orders in columns
    #    where m = 0,1,...l
    # colat: 1D array of the colatitudes at each point of field X
    #        along vertical dimension (in radians) 
    # lons: longitudes of grid (in radians from 0 to 2pi)
    # lmax: the maximum degree of spherical harmonic to convert
    # dphi: rotate field by dphi (radians) across longitudinal direction

    # outputs:
    # X: 2D array of the global field rotated
    # Ulm_rot: rotated Ulms that will map back onto 0 - 360

    import numpy as np
    import scipy.special as sp
    import math as ms

    MM=len(colat)
    NN=len(lons)
    X=np.zeros((MM,NN))
    Y=np.zeros((MM,NN),dtype=complex)
    Ulm_rot=np.copy(U)
    

    for l in range(lmax+1):
        # first the 0-th order
        for i in range(MM):
            #Y[i,:]=sp.sph_harm(0,l,lons,colat[i])
            m=0
            FACT=np.sqrt((2.*l+1.)*ms.factorial(l-m)/ \
                                        ms.factorial(l+m))            
	    #FACT=np.sqrt(2.)
            Y[i,:]=sp.lpmv(m,l,np.cos(colat[i]))*FACT* \
                np.exp(1.j*m*lons)
        X[:,:]=X[:,:]+np.real(U[l,0])*np.real(Y)

    for l in range(1,lmax+1):
        for m in range(1,l+1):
            for i in range(MM):
            # form spherical harmonics
                FACT= np.sqrt((2.*l+1.)*ms.factorial(l-m)/ \
                                        ms.factorial(l+m))            
                #FACT= np.sqrt(2.)
		Y[i,:]=sp.lpmv(m,l,np.cos(colat[i]))*FACT* \
                    np.exp(1.j*m*lons)

#                Y[i,:]=sp.sph_harm(m,l,lons,colat[i])
            # perform summation    
            X[:,:]=X[:,:]+2.*np.real(U[l,m]*Y*\
                                         np.exp(1.j*m*dphi)) # rotate   
            Ulm_rot[l,m]=U[l,m]*np.exp(1.j*m*dphi)
                
    return X,Ulm_rot


def norm_sph(colat,lon,l,m):
    # by Harriet Lau; 2 Jul 2013
    # produces spherical harmonic in (theta,phi) coordinates following the
    # the normalization given in J. X. Mitrovica's Sea Level equation
    # This is *not* the same normalization as in seismology literature.
    # {Sea level Ylm} = (4pi)**0.5{seismology Ylm}

    # inputs:
    # colat: 1D array/scalar of colatitude nodes in which to calculate
    #        spherical harmonics on (0 - pi radians)
    # lon: 1D array/scalar of longitude nodes in which to calculate 
    #      spherical harmonics
    # l: degree to calculate
    # m: order to calculate

    # outputs:
    # Y: 2D array of spherical harmonics at colat and lon points

    import numpy as np
    import scipy.special as sp
    import math as ms

    MM=len(colat)
    NN=len(lon)
    Y=np.zeros((MM,NN),dtype=complex)

    FACT=(-1.)**m * np.sqrt((2.*l+1.)*ms.factorial(l-m)/ \
                                        ms.factorial(l+m))            
    #FACT= np.sqrt(2.)
    for i in range(MM):
        Y[i,:]=sp.lpmv(m,l,np.cos(colat[i]))*FACT* \
            np.exp(1.j*np.float(m)*lon)

    return Y

def norm_sph_sites(colat,lon,l,m):
    # by Harriet Lau; 2 Jul 2013
    # produces spherical harmonic in (theta,phi) coordinates following the
    # the normalization given in J. X. Mitrovica's Sea Level equation
    # This is *not* the same normalization as in seismology literature.
    # {Sea level Ylm} = (4pi)**0.5{seismology Ylm}

    # inputs:
    # colat: 1D array/scalar of colatitude nodes in which to calculate
    #        spherical harmonics on (0 - pi radians) - of N length
    # lon: 1D array/scalar of longitude nodes in which to calculate 
    #      spherical harmonics - of N length 
    # l: degree to calculate
    # m: order to calculate

    # outputs:
    # Y: 1D list of Ylms at N sites for given colat/lon pairs

    import numpy as np
    import scipy.special as sp
    import math as ms

    MM=len(colat)
    Y=np.zeros(MM,dtype=complex)

    FACT=(-1.)**m * np.sqrt((2.*l+1.)*ms.factorial(l-m)/ \
                                        ms.factorial(l+m))            
    #FACT= np.sqrt(2.)
    for i in range(MM):
        Y[i]=sp.lpmv(m,l,np.cos(colat[i]))*FACT* \
            np.exp(1.j*np.float(m)*lon[i])

    return Y

def real_sph(colat,lon,l,m):
    # by Harriet Lau; 25 Sep 2014
    # produces real spherical harmonic in (theta,phi) coordinates following the
    # the normalization given in
    # http://www.boost.org/doc/libs/1_49_0/libs/math/ ...
    # doc/sf_and_dist/html/math_toolkit/special/sf_poly/sph_harm.html

    # inputs:
    # colat: 1D array/scalar of colatitude nodes in which to calculate
    #        spherical harmonics on (0 - pi radians)
    # lon: 1D array/scalar of longitude nodes in which to calculate 
    #      spherical harmonics
    # l: degree to calculate
    # m: order to calculate

    # outputs:
    # Y: 2D array of spherical harmonics at colat and lon points

    import numpy as np
    import scipy.special as sp
    import math as ms

    MM=len(colat)
    NN=len(lon)
    Y=np.zeros((MM,NN))

    mm=np.abs(m)
    FACT=(-1.)**mm * np.sqrt((2.*l+1.)/(4.*np.pi)*ms.factorial(l-mm)/ \
                                        ms.factorial(l+mm))            


    

    if (m<0):
        prefac=np.cos(m*lon)#*np.sqrt(2.) HL
    elif (m==0):
        prefac=1.
    else:
        prefac=np.sin(m*lon)#*np.sqrt(2.)

    for i in range(MM):
        Y[i,:]=sp.lpmv(mm,l,np.cos(colat[i]))*FACT*\
            prefac

    return Y

def grad_sph(colat,lon,l,m,r):
    # by Harriet Lau; Oct 2014 
    # calculates the lateral gradient of a spherical harmonic of choice
    # inputs:
    # colat: 1D array/scalar of colatitude nodes (0 - pi radians)
    # lon: 1D array/scalar of longitude nodes in which to calculate
    #      spherical harmonics
    # l: degree to calculate
    # m: order to calculate
    # r: radius
    # outputs:
    # grad_Ylm1: array/scalar of grad[Ylm] in colat direction 
    #            (at points specified by colat and lon)
    # grad_Ylm2: array/scalar of grad[Ylm] in long direction
    #            (at points specified by colat and lon)
    # Using normalization of JXM.
    

    import numpy as np
    import scipy.special as sp
    import math as ms
    
    if (l==0):
        print 'GRAD_SPH CANNOT TAKE L=0, ZERO GRADIENT'

    MM=len(colat)
    NN=len(lon)
    grad_Ylm1=np.zeros((MM,NN),dtype=complex)
    grad_Ylm2=np.zeros((MM,NN),dtype=complex)

    FACT=(-1.)**m * np.sqrt((2.*l+1.)*ms.factorial(l-m)/ \
                                        ms.factorial(l+m))            
    FACT1a=(-1.)**(m+1) * np.sqrt((2.*l+1.)*ms.factorial(l-(m+1))/ \
                                        ms.factorial(l+(m+1)))            
    FACT1b=0.5*np.sqrt((l-m)*(l+m+1))/r
    FACT2a=(-1.)**(m-1) * np.sqrt((2.*l+1.)*ms.factorial(l-(m-1))/ \
                                        ms.factorial(l+(m-1)))            
    FACT2b=0.5*np.sqrt((l+m)*(l-m+1))/r

    for i in range(MM):
        grad_Ylm1[i,:]=(FACT1a*FACT1b*sp.lpmv(m+1,l,np.cos(colat[i]))-\
                            FACT2a*FACT2b*sp.lpmv(m-1,l,np.cos(colat[i])))*\
                            np.exp(1.j*m*lon)
        grad_Ylm2[i,:]=FACT*sp.lpmv(m,l,np.cos(colat[i]))*\
            1.j*m*np.exp(1.j*m*lon)/(r*np.sin(colat[i]))
        
    return grad_Ylm1,grad_Ylm2
