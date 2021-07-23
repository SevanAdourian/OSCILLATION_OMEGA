def list_to_grid(X,n):
    # by Harriet Lau; 6 Aug 2013
    # converts output X from Latychev's finite volume GIA code
    # X is a list of values that vary in both latitude and 
    # longitude.  For latitude this is a Gauss Legendre grid.

    # X is converted to a grid suitable for contour plots
    
    # inputs:
    # X: 1D array of values of a global field
    # n: degree of Gauss-Legendre nodes
    # length of X must be 2*n**2

    # outputs:
    # colat: colatitudes of grid nodes (radians); n-long 
    # lons: longitudes of grid nodes (radians); 2n-long
    # XX: 2D grid (n x 2n) 

    import scipy.special as sp
    import numpy as np
    import numpy.polynomial.legendre as leg

    [colats,dx]=leg.leggauss(n)
    colats=np.arccos(colats)
    lons=np.linspace(0,2.*np.pi,2*n+1)
    lons=lons[0:-1] # so we do not repeat the same longitude (0,2pi)
    XX=np.reshape(X,(n,2*n))

    return dx,colats,lons,XX
