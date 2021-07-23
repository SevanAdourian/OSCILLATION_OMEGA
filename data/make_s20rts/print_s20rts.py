#! /usr/bin/env python

import numpy as np
import s24rts as s24

# import depths from PREM:
r = np.loadtxt('isoprem808.md',skiprows=3,usecols=(0),unpack=True)
icmb=330
r = r[icmb:len(r)]/1.e3
ndepth = len(r)

rout = np.zeros((ndepth,2))
rout[:,0] = range(ndepth)
rout[:,1] = r
np.savetxt('rad.dat',rout,fmt='%f',header='lay number; radius (km)')


for i in range(ndepth):
   print 'layer',i+1,'of',ndepth
   d = 6371.-r[i]
   [ulms,junk,junk] = s24.s24rts('S20',d,0,0,1,0)
   ulms = 0.4*ulms*np.sqrt(4.*np.pi) # right normalization for seismology
   # 0.4 scaling to density

   rout = np.real(ulms)
   iout = np.imag(ulms)
   fre = "rho_ulm/rho_ulm_re_lay%03d.dat" % i
   fim = "rho_ulm/rho_ulm_im_lay%03d.dat" % i
   # fim = 'coucou'
   print fim
   np.savetxt(fre,rout,fmt='%e')
   np.savetxt(fim,iout,fmt='%e')
   # np.savetxt('rho_ulm_im_lay'+str(i)+'.dat',iout,fmt='%e')


 





