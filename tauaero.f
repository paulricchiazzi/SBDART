c file:                  tauaero.f
c
c
c external routines:     tauaero,relhum
c
c required routines:     locate, numset
c
c internal routines:     aeroden,aerzstd,phaerw,aerbwi,aestrat,
c                        usraer,stdaer, aervint,aeread

c
c=======================================================================
      module aeroblk

      use params, only: mxly, maxmom, zip, kr
      implicit none

      integer, parameter ::
     &     naerz=5,             ! number of stratospheric aerosol layers
     &     naerb=150,           ! max number of BLA wavelengths
     &     naerw=47,            ! number of wavelengths in standard models
     &     naerwm=naerb*maxmom  ! max pmaer buffer size


      ! quantities set in main program:

      real(kr) ::
     &     zaer(naerz)=0. ,     ! stratospheric aerosol layer altitudes
     &     taerst(naerz)=0.,    ! stratospheric aerosol layer optical depth
     &     zbaer(mxly)=zip,     ! user definded layer altitudes (BLA only)
     &     dbaer(mxly)=zip,     ! user defined density profile (BLA only)
     &     vis=zip,             ! visibility at 0.55um (BLA only)
     &     tbaer=zip,           ! BLA optical depth at 0.55um
     &     abaer=0.,            ! BLA lambda exponent 
     &     wlbaer(naerb)=zip,   ! BLA wavelengths 
     &     qbaer(naerb)=zip,    ! BLA extinction efficiency
     &     wbaer(naerb)=zip,    ! BLA single scattering albedo
     &     gbaer(naerb)=zip,    ! BLA asymmetry factor
     &     pmaer(naerwm)=zip,   ! BLA phase function moments
     &     rhaer=zip            ! relative humidity (BLA only)

      integer ::
     &     iaer=0,              ! BLA aerosol type 
     &     imoma=3,             ! flag to control use of pmaer (BLA only)
     &     jaer(naerz)=0        ! stratospheric aerosol types

      ! the following quantities not accessed by main program

      integer ::
     &     nzbaer,              ! number of user specified BLA levels
     &     ndbaer,              ! number of user specified BLA densities
     &     nwlbaer=naerw,       ! number of user specified BLA wavelengths
     &     npmaer               ! number of user specified BLA pmaer moments

      real(kr), dimension(naerb) ::
     &     aerext(naerb),       ! BLA extinction
     &     aerabs(naerb),       ! BLA absorption
     &     aerasm(naerb)        ! BLA asymmetry

      real(kr), dimension(naerw) :: ! wavelengths for standard BLA and strat models
     &     awl=(/
     &   .2000,   .3000,   .3371,   .5500,   .6943,  1.0600,  1.5360,
     &  2.0000,  2.2500,  2.5000,  2.7000,  3.0000,  3.3923,  3.7500,
     &  4.5000,  5.0000,  5.5000,  6.0000,  6.2000,  6.5000,  7.2000,
     &  7.9000,  8.2000,  8.7000,  9.0000,  9.2000, 10.0000, 10.5910,
     & 11.0000, 11.5000, 12.5000, 14.8000, 15.0000, 16.4000, 17.2000,
     & 18.5000, 21.3000, 25.0000, 30.0000, 40.0000, 50.0000, 60.0000,
     & 80.0000, 100.000, 150.000, 200.000, 300.000/)
                  
      contains
c-----------------------------------------------------------------------

      subroutine aerzstd
c
c purpose 
c       copy standard vertical profile into zbaer and dbaer and change
c       nzbaer to number of elements in standard profile
c
c the standard vertical distributions of the boundary layer aerosol 
c density is based on the 5s vertical profile models for 5 and 23 km 
c visibility.
c
c above 5 km, the aden05 and aden23 models are the same 
c below 5 km, the models differ as follows;
c aden05     0.99 km scale height (94% of extinction occurs below 5 km)
c aden23     1.45 km scale height (80% of extinction occurs below 5 km)
c    

      use params, only: one, zero

      implicit none

      integer, parameter ::  mz=33 ! number of standard BLA layer altitudes
      real(kr), dimension(mz) :: alt=(/
     &        0.0,        1.0,        2.0,        3.0,        4.0,  
     &        5.0,        6.0,        7.0,        8.0,        9.0,
     &       10.0,       11.0,       12.0,       13.0,       14.0, 
     &       15.0,       16.0,       17.0,       18.0,       19.0,
     &       20.0,       21.0,       22.0,       23.0,       24.0, 
     &       25.0,       30.0,       35.0,       40.0,       45.0,
     &       50.0,       70.0,      100.0/)
        real(kr), dimension(mz) :: aden05=(/
     &  1.378E+04,  5.030E+03,  1.844E+03,  6.731E+02,  2.453E+02,      
     &  8.987E+01,  6.337E+01,  5.890E+01,  6.069E+01,  5.818E+01,
     &  5.675E+01,  5.317E+01,  5.585E+01,  5.156E+01,  5.048E+01,
     &  4.744E+01,  4.511E+01,  4.458E+01,  4.314E+01,  3.634E+01,
     &  2.667E+01,  1.933E+01,  1.455E+01,  1.113E+01,  8.826E+00,
     &  7.429E+00,  2.238E+00,  5.890E-01,  1.550E-01,  4.082E-02,
     &  1.078E-02,  5.550E-05,  1.969E-08/)
        real(kr), dimension(mz) :: aden23=(/
     &  2.828E+03,  1.244E+03,  5.371E+02,  2.256E+02,  1.192E+02,
     &  8.987E+01,  6.337E+01,  5.890E+01,  6.069E+01,  5.818E+01,
     &  5.675E+01,  5.317E+01,  5.585E+01,  5.156E+01,  5.048E+01,
     &  4.744E+01,  4.511E+01,  4.458E+01,  4.314E+01,  3.634E+01,
     &  2.667E+01,  1.933E+01,  1.455E+01,  1.113E+01,  8.826E+00,
     &  7.429E+00,  2.238E+00,  5.890E-01,  1.550E-01,  4.082E-02,
     &  1.078E-02,  5.550E-05,  1.969E-08/)

      integer :: i
      real(kr) :: wtv

      if(mz.gt.mxly) then
        print *,'mz gt mxly in aerzstd'
        stop
      endif

      nzbaer=mz

      if(vis.le.0.) then
        wtv=(1./vis-1/5.)/(1./23.-1./5.)
      else
        wtv=1.                  ! default to 23 km vertical profile
      endif
      wtv=max(zero,min(one,wtv))

      do i=1,mz
        zbaer(i)=alt(i)
        dbaer(i)=aden05(i)*(1.-wtv)+aden23(i)*wtv
      enddo
      return
      end subroutine aerzstd
      
c=======================================================================
      subroutine phaerw(w,p)
c
c purpose:
c   Interpolate phase function moments to wavelength w.
c   Must use linear interpolation because moments can be pos or neg.
c   if w>wmax then p(w)=ptable(wmax)
c   if w<wmin then p(w)=ptable(wmin)
c
c input
c   w      wavelength in microns
c
c output
c   p      legendre moments of scattering phase function (0:nstrms)
c          due to boundary layer aerosol
c
      use params, only: one, zero,kr
      implicit none
      integer :: l,j,i
      real(kr) :: wt,w,p(*)

      call locate(wlbaer,nwlbaer,w,l)
      wt=(w-wlbaer(l))/(wlbaer(l+1)-wlbaer(l))
      wt=max(zero,min(wt,one))     ! no extrapolation

      do i=1,npmaer
        j=l+(i-1)*nwlbaer
        p(i)=pmaer(j)*(1.-wt)+pmaer(j+1)*wt
      enddo

      return
      end subroutine phaerw

c=======================================================================
      subroutine aerbwi(wl,extinc,wa,ga)
c 
c purpose
c     At a given wavelength find the extinction efficiency,
c     single scattering albedo and asymmetry factor of boundary
c     layer aerosols by interpolation
c
c input
c  wl         wavelength
c 
c output
c  extinc     extinction efficiency
c  wa         single scattering albedo
c  ga         asymmetry factor

      use params, only: one, zero, kr
      implicit none

      real(kr),intent(in) :: wl
      real(kr),intent(out) :: extinc,wa,ga
      real(kr) :: wt,absorp

      integer :: l
      character*(80) errmes(10)
c
      LOGICAL, parameter :: debug=.false.


c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      IF (debug) write(*,*) 'db_aerbwi.1',wl

      if(iaer.eq.0) then   ! no boundary layer aerosol
        extinc=0.
        wa=0.
        ga=0.
        return
      endif        

c     print '(a/(10es11.3))','wlbaer',wlbaer(1:nwlbaer)
c     print '(a/(10es11.3))','aerext',aerext(1:nwlbaer)
c     stop !!!!!!!!!!!!!!!!!!!!

      call locate(wlbaer,nwlbaer,wl,l)

      if(wl.le.wlbaer(1)) then 
        extinc=aerext(1)*(wlbaer(1)/wl)**abaer
        wa=1.-(aerabs(1)/aerext(1))
        ga=aerasm(1)
      elseif(wl.ge.wlbaer(nwlbaer)) then 
        extinc=aerext(nwlbaer)*(wlbaer(nwlbaer)/wl)**abaer
        wa=1.-(aerabs(nwlbaer)/aerext(nwlbaer))
        ga=aerasm(nwlbaer)
      else
        wt=log(wl/wlbaer(l))/log(wlbaer(l+1)/wlbaer(l))

        extinc=aerext(l)*(aerext(l+1)/aerext(l))**wt
        if(aerabs(l).gt.0. .and. aerabs(l+1).gt.0. ) then 
          absorp=aerabs(l)*(aerabs(l+1)/aerabs(l))**wt
        else
          absorp=aerabs(l)*(1.-wt)+aerabs(l+1)*wt
        endif
        if(extinc.gt.0.) wa=max(zero,min(one-absorp/extinc,one))
        ga=(1.-wt)*aerasm(l)+wt*aerasm(l+1)  
        ! linear interpolation because aerasm can be negative
      endif
      !write(*,*) 'wlbaer: ',wlbaer(1),wlbaer(2)
      !write(*,*) 'extinc: ',extinc
      !write(*,*) 'aerext: ',aerext(1),aerext(2)
      !write(*,*) 'wl,l ',wl,l

      return
      end subroutine aerbwi

c=======================================================================

      subroutine aestrat(ja,wl,qa,wa,ga)

c 
c purpose
c     compute extinction efficiency, single scattering albedo and
c     asymmetry factor of STRATOSPHERIC aerosols by interpolation
c
c input
c  ja         stratospheric aerosol type
c  wl         wavelength
c 
c output
c  qa         extinction efficiency of stratospheric aerosol of type ja
c  wa         single scattering albedo
c  ga         asymmetry factor
c

      use params, only: zero, one, kr
      implicit none
      real(kr) :: wl,wt,qa,wa,ga,absorp
      integer :: l,ja
      real(kr), dimension(naerw,3,4) :: aerstr

                                ! background extinction
      data aerstr(1:naerw,1,1) /       
     & 2.0752e0,  1.8656e0,  1.7246e0,  1.0000e0,  7.0156e-1, 3.0170e-1,
     & 1.1440e-1, 5.1225e-2, 3.4285e-2, 2.3475e-2, 1.6878e-2, 6.6506e-2,
     & 1.0943e-1, 8.9653e-2, 6.7609e-2, 5.2855e-2, 6.7496e-2, 5.7975e-2,
     & 4.2471e-2, 2.4176e-2, 4.6102e-2, 1.2339e-1, 1.7699e-1, 1.2389e-1,
     & 9.0220e-2, 8.5793e-2, 3.2838e-2, 2.6528e-2, 5.0703e-2, 1.9471e-2,
     & 1.1710e-2, 1.6106e-2, 1.7716e-2, 3.9533e-2, 3.7954e-2, 5.4871e-3,
     & 8.8409e-3, 1.2289e-3, 1.0647e-3, 3.3151e-3, 4.5164e-3, 4.1496e-3,
     & 3.2801e-3, 2.4481e-3, 1.0714e-3, 5.0381e-4, 1.2101e-4/
                                ! background absorption
      data aerstr(1:naerw,2,1) /
     & 2.4347e-7, 1.4949e-7, 1.3020e-7, 6.8716e-8, 1.0053e-7, 4.2384e-6,
     & 2.3075e-4, 1.4889e-3, 1.8377e-3, 3.3645e-3, 4.6858e-3, 5.9424e-2,
     & 1.0129e-1, 8.2396e-2, 6.3621e-2, 5.0327e-2, 6.5798e-2, 5.6149e-2,
     & 4.0939e-2, 2.3226e-2, 4.5800e-2, 1.2276e-1, 1.7568e-1, 1.2216e-1,
     & 8.9058e-2, 8.4861e-2, 3.1979e-2, 2.6097e-2, 5.0244e-2, 1.9042e-2,
     & 1.1497e-2, 1.6024e-2, 1.7639e-2, 3.9452e-2, 3.7847e-2, 5.4251e-3,
     & 8.8160e-3, 1.2118e-3, 1.0579e-3, 3.3122e-3, 4.5150e-3, 4.1489e-3,
     & 3.2799e-3, 2.4481e-3, 1.0713e-3, 5.0381e-4, 1.2101e-4/
                                ! background asymmetry
      data aerstr(1:naerw,3,1) /
     & 0.6749, 0.6943, 0.6991, 0.6846, 0.6572, 0.5861, 0.4965, 0.4164,
     & 0.3772, 0.3385, 0.3069, 0.2599, 0.2234, 0.2028, 0.1554, 0.1291,
     & 0.1057, 0.0962, 0.0909, 0.0802, 0.0595, 0.0458, 0.0413, 0.0479,
     & 0.0483, 0.0451, 0.0504, 0.0379, 0.0346, 0.0365, 0.0273, 0.0180,
     & 0.0174, 0.0144, 0.0151, 0.0141, 0.0094, 0.0075, 0.0049, 0.0031,
     & 0.0020, 0.0014, 0.0008, 0.0005, 0.0002, 0.0001, 0.0001/
                                ! aged volcanic extinction
      data aerstr(1:naerw,1,2) /
     & 1.14880, 1.19171, 1.18013, 1.00, .84873, .53019, .27968, .14551,
     &  .11070, .08633, .07184, .06076, .04506, .03399, .02095, .01538,
     &  .01266, .01019, .00994, .01044, .01361, .01791, .02278, .02918,
     &  .03108, .03234, .03456, .03184, .02772, .02475, .01715, .01563,
     &  .01665, .01646, .01734, .01772, .01076, .01051, .01133, .01329,
     &  .01492, .01277, .00766, .00562, .00318, .00231, 0.0/
                                ! aged volcanic absorption
      data aerstr(1:naerw,2,2) /
     &  .44816, .11259, .08500, .05272, .04082, .02449, .01487, .01019,
     &  .00867, .00842, .00842, .00949, .00741, .00487, .00316, .00335,
     &  .00399, .00449, .00525, .00665, .01114, .01652, .02177, .02437,
     &  .02506, .02658, .03006, .02861, .02513, .02285, .01620, .01532,
     &  .01633, .01620, .01709, .01741, .01057, .01038, .01127, .01329,
     &  .01492, .01277, .00766, .00562, .00318, .00231, 0.0/
                                ! aged volcanic asymmetry
      data aerstr(1:naerw,3,2) /
     &  .8272, .7148, .7076, .6978, .6886, .6559, .6062, .5561, .5255,
     &  .4958, .4729, .4401, .4015, .3699, .3125, .2773, .2472, .2173,
     &  .2054, .1908, .1623, .1348, .1233, .1615, .1757, .1712, .1521,
     &  .1326, .1230, .1081, .0801, .0528, .0514, .0461, .0446, .0449,
     &  .0415, .0330, .0198, .0097, .0044, .0032, .0020, .0013, .0006,
     &  .0004, .0000/
                                ! fresh volcanic extinction
      data aerstr(1:naerw,1,3) /
     &  .88715, .92532, .94013, 1.0000, 1.03013, 1.05975, 1.01171,
     &  .88677, .82538, .76361, .71563, .67424, .60589, .55057, .45222,
     &  .37646, .32316, .25519, .22728, .20525, .17810, .14481, .14152,
     &  .37639, .44551, .44405, .42222, .36462, .32551, .27519, .16728,
     &  .10627, .10861, .10886, .11665, .13127, .10108, .08557, .06411,
     &  .05741, .05531, .04707, .02792, .02028, .01136, .00820, 0.0/
                                ! fresh volcanic absorption
      data aerstr(1:naerw,2,3) /
     &  .41582, .22892, .19108, .14468, .12475, .09158, .06601, .04943,
     &  .04367, .04342, .04399, .05076, .04133, .02829, .01924, .01981,
     &  .02297, .02475, .02778, .03411, .05335, .07133, .08816, .15342,
     &  .18506, .19354, .20791, .18449, .16101, .13759, .08456, .06886,
     &  .07278, .07367, .07956, .08785, .06032, .05747, .05133, .05323,
     &  .05453, .04657, .02773, .02020, .01135, .00820, 0.0/
                                ! fresh volcanic asymmetry
      data aerstr(1:naerw,3,3) /
     &  .9295, .8115, .7897, .7473, .7314, .7132, .7113, .7238, .7199,
     &  .7165, .7134, .6989, .6840, .6687, .6409, .6325, .6199, .6148,
     &  .6142, .6072, .5853, .5632, .5486, .4753, .4398, .4329, .4091,
     &  .4105, .4120, .4136, .4140, .3637, .3577, .3344, .3220, .3052,
     &  .2957, .2564, .2055, .1229, .0632, .0483, .0321, .0216, .0103,
     &  .0059, .0000/
                                ! meteoric extinction
      data aerstr(1:naerw,1,4) /

     & 1.05019, 1.05880, 1.05259, 1.00, .94949, .81456, .66051, .54380,
     &  .49133, .44677, .41671, .38063, .34778, .32804, .29722, .27506,
     &  .25082, .22620, .21652, .20253, .17266, .14905, .14234, .14082,
     &  .15057, .16399, .23608, .24481, .27791, .25076, .15272, .09601,
     &  .09456, .14576, .12373, .18348, .12190, .12924, .08538, .04108,
     &  .04714, .04069, .02480, .01789, .00980, .00693, 0.0/
                                ! meteoric absorption
      data aerstr(1:naerw,2,4) /
     &  .00063, .00152, .00184, .00506, .00791, .01829, .03728, .06158,
     &  .07538, .08943, .10051, .11614, .13310, .14348, .14633, .13728,
     &  .12462, .11184, .10709, .10076, .09006, .08734, .09000, .10304,
     &  .11905, .13437, .19551, .20095, .22494, .18418, .09285, .06665,
     &  .06823, .12329, .10551, .16184, .09835, .10582, .06759, .03247,
     &  .04405, .03816, .02327, .01696, .00946, .00677, 0.0/
                                ! meteoric asymmetry
      data aerstr(1:naerw,3,4) /
     &  .7173, .7039, .7020, .6908, .6872, .6848, .6891, .6989, .7046,
     &  .7099, .7133, .7159, .7134, .7058, .6827, .6687, .6583, .6513,
     &  .6494, .6475, .6467, .6496, .6506, .6461, .6334, .6177, .5327,
     &  .5065, .4632, .4518, .5121, .5450, .5467, .4712, .4853, .3984,
     &  .4070, .3319, .3427, .3766, .3288, .2969, .2808, .2661, .2409,
     &  .2098, .0000/

      call locate(awl,naerw,wl,l)
      if(wl.le.awl(1)) then 
        qa=aerstr(1,1,ja)*(awl(1)/wl)**abaer
        wa=1.-(aerstr(1,2,ja)/aerstr(1,1,ja))
        ga=aerstr(l,3,ja)
      elseif(wl.ge.awl(naerw)) then 
        qa=aerstr(naerw,1,ja)*(awl(1)/wl)**abaer
        wa=1.-(aerstr(naerw,2,ja)/aerstr(naerw,1,ja))
        ga=aerstr(naerw,3,ja)
      else
        wt=log(wl/awl(l))/log(awl(l+1)/awl(l))
        qa=aerstr(l,1,ja)*(aerstr(l+1,1,ja)/aerstr(l,1,ja))**wt
        absorp=aerstr(l,2,ja)*(aerstr(l+1,2,ja)/aerstr(l,2,ja))**wt
        if(qa.gt.0.) wa=max(zero,min(one-absorp/qa,one))
        ga=(1.-wt)*aerstr(l,3,ja)+wt*aerstr(l+1,3,ja)  
        ! linear interpolation because aerasm can be negative
      endif
      !write(*,*) 'awl: ',awl(1),awl(2)
      !write(*,*) 'qa: ',qa
      !write(*,*) 'aerext: ',aerext(1),aerext(2)
      !write(*,*) 'wl,l ',wl,l

      return
      end subroutine aestrat

c=======================================================================

      subroutine usraer(q55)
c
c purpose 
c       Replaces data in wlbaer,aerext,aerabs,aerasm with user specified
c       aerosol spectral data.  
c
c       A similar function is performed by subroutine STDAER, which
c       fills the same variables with spectral properties of standard
c       aerosol models
c
c input:
c
c  wlbaer    wavelengths at which aerosol model is specified
c  qbaer    aerosol extinction efficiency
c  wbaer    aerosol single scattering albedo
c  gbaer    aerosol asymmetry factor
c  abaer    power law slope used to extrapolate optical depth
c           outside of wlbaer range  tauaer=(wlbaer(j)/wl)**abaer
c           where j is 1 or n (n = total number of assigned values
c           of wlbaer). 
c  pmaer    user defined aerosol phase function (may be undefined)
c
c input/output:
c
c  q55      extinction efficiency at 0.55 um
c
c output:
c
c  imoma    >0 => aerosol phase function from subroutine getmom
c            0 => user defined aerosol phase function
c

c check input parameters

      use params, only: wl55, zero, kr

      implicit none
      integer, external :: numset
      integer :: nqbaer, nwbaer, ngbaer, ne, i, j
      real(kr) :: f,q55
      character (len=80) :: errmes(10)

      nwlbaer=numset(zip,wlbaer,naerb)
      nqbaer=numset(zip,qbaer,naerb)
      nwbaer=numset(zip,wbaer,naerb)
      ngbaer=numset(zip,gbaer,naerb)
      npmaer=numset(zip,pmaer,naerwm)

      ne=0

      if(nwlbaer.eq.0) then         ! wavelength points not set

        qbaer(1)=1.
        nqbaer=1
        if(nwbaer.ne.1) then
          ne=ne+1
          errmes(ne)='specify one value of wbaer when wlbaer not set'
        endif

      elseif(nwlbaer.eq.1) then     ! only one wavelength point

        if(nqbaer.gt.1) then 
          ne=ne+1
          errmes(ne)='number of elements must match: wlbaer, qbaer'
        elseif(nqbaer.eq.0) then
          qbaer(1)=1.
        endif          
        if(nwbaer.ne.1) then
          ne=ne+1
          errmes(ne)='number of elements must match: wlbaer, wbaer'
        endif

      else                          ! many wavelength points

        if(nwlbaer.ne.nqbaer) then
          ne=ne+1
          errmes(ne)='number of elements must match: wlbaer, qbaer'
        endif
        if(nwlbaer.ne.nwbaer) then
          ne=ne+1
          errmes(ne)='number of elements must match: wlbaer, wbaer'
        endif
        if(ngbaer.eq.0) then
          if(npmaer.ge.1) then 
            if(mod(npmaer,nwlbaer).ne.0) then
              ne=ne+1
              errmes(ne)='incorrect number of phase function moments'
            else
              npmaer=npmaer/nwlbaer
            endif
          else
            ne=ne+1
            errmes(ne)='must specify either gbaer or pmaer'
          endif
        else
          if(ngbaer.ne.nwbaer) then
            ne=ne+1
            errmes(ne)='number of elements must match: wlbaer, gbaer'
          endif
        endif
      
      endif

      if((ngbaer.eq.0).eqv.(npmaer.eq.0).and.imoma.eq.3) then 
        ne=ne+1
        errmes(ne)='must specify either gbaer or pmaer, not both'
      endif

      if(ne.gt.0) then 
        write(*,*) 'Error in user specified aerosols (iaer=5)'
        write(*,'(/,1x,5a8)')
     &          'nwlbaer','nqbaer','nwbaer','ngbaer','npmaer'
        write(*,'(5i8,/)')  nwlbaer,nqbaer,nwbaer,ngbaer,npmaer
        do i=1,ne
          write(*,'(2a)') 'Error in USRAER -- ',errmes(i)
        enddo
        stop
      endif
      
c duplicate input to two wavelength points if a 
c single wavelength specified in the input file

      !write(*,*) 'db usraer 1'
      if(nwbaer.eq.1) then  
        if(nwlbaer.eq.0) wlbaer(1)=wl55
        wlbaer(2)=2*wlbaer(1)
        nwlbaer=2
        qbaer(2)=qbaer(1)*(wlbaer(1)/wlbaer(2))**abaer
        wbaer(2)=wbaer(1)
        gbaer(2)=gbaer(1)

c duplicate phase function entries

        !write(*,*) 'db usraer 2'
        if(npmaer.gt.0) then
          pmaer(2:2*npmaer:2)=pmaer(1:npmaer)
          pmaer(1:2*npmaer-1:2)=pmaer(2:2*npmaer:2)
        endif
      endif

      !write(*,'(a10,/,10es11.3)')  'wlbaer',(wlbaer(i),i=1,2)
      !write(*,'(a10,/,10es11.3)')  'qbaer',(qbaer(i),i=1,2)
        
      !write(*,*) 'db usraer 3'
      if(nwlbaer.gt.naerb) then
        write(*,*) 'Too many spectral points in usraer'
        stop
      endif
      do i=1,nwlbaer 
        wlbaer(i)=wlbaer(i)
        aerext(i)=qbaer(i)
        aerabs(i)=(1.-wbaer(i))*aerext(i)
        aerasm(i)=gbaer(i)
      enddo

      !write(*,*) 'db usraer 4'
      !write(*,'(a10/10es11.3)')  'wlbaer',(wlbaer(i),i=1,nwlbaer)
      !write(*,'(a10/10es11.3)')  'aerext',(aerext(i),i=1,nwlbaer)
      !write(*,'(a10/10es11.3)')  'aerabs',(aerabs(i),i=1,nwlbaer)
      !write(*,'(a10/10es11.3)')  'aerasm',(aerasm(i),i=1,nwlbaer)
      !write(*,*) 'pmaer'
      !do j=1,npmaer
      !  write(*,'(10es11.3)') (pmaer((j-1)*nwlbaer+i),i=1,nwlbaer)
      !enddo
        
      call locate(wlbaer,nwlbaer,wl55,j)

      !write(*,*) 'db usraer 5'
      !write(*,*) 'wlbaer ',j,wlbaer(j),wlbaer(j+1)
      f=log(wl55/wlbaer(j))/log(wlbaer(j+1)/wlbaer(j))
      !write(*,*) 'f ',f
      if(wl55.lt.wlbaer(1)) then
        !write(*,*) 'db usraer 5.1',abaer
        q55=qbaer(1)*(wlbaer(1)/wl55)**abaer
      elseif(wl55.gt.wlbaer(nwlbaer)) then
        !write(*,*) 'db usraer 5.2',abaer
        q55=qbaer(nwlbaer)*(wlbaer(nwlbaer)/wl55)**abaer
      else
        !write(*,*) 'db usraer 5.3',qbaer(j),f
        q55=qbaer(j)*(qbaer(j+1)/qbaer(j))**f
      endif

c if defined, save aerosol phase function for use in phaerw

      if(npmaer.gt.0) imoma=0

      return

      end subroutine usraer

c===========================================================================
      subroutine stdaer(humid)
c
c purpose 
c   Replaces data in wlbaer,aerext,aerabs,aerasm with standard
c   extinction, absorption, and asymmetry parameters for boundary
c   layer aerosols.  The rural, urban, oceanic or tropospheric
c   aerosols models are interpolated on relative humidity.  
c
c   A similar function is performed by entry USRAER, which fills the
c   same variables with spectral properties of user defined aerosol model

c
c input:
c
c   humid    relative humidity at the surface (0.0 - 1.0)
c
      
      use params, only: wl55, zero, one, kr
      implicit none

      integer :: j,i

      real(kr) :: rhzone(4)=(/0.,.7,.8,.99/),
     &        tiny=.00000001,
     &        rhum,humid,wt,ex1=0.0_kr,ex2=0.0_kr,
     &        ab1=0.0_kr,ab2=0.0_kr,g1=0.0_kr,g2=0.0_kr

      real(kr), dimension(naerw,4) ::
     &     rure,rura,rurg,      ! rural (four humidity values)
     &     urbe,urba,urbg,      ! urban
     &     ocne,ocna,ocng,      ! oceanic
     &     troe,troa,trog       ! tropospheric

c
c
c   altitude regions for aerosol extinction coefficients
c
c         aerosol extinction normalized to one at 0.55 microns
c
c         0-2km
c           rure,rura,rurg=rural extinction, absorption, asymmetry
c           urbe,urba,urbg=urban extinction, absorption, asymmetry
c           ocne,ocna,ocng=maritime extinction, absorption, asymmetry
c           troe,troa,trog=tropspheric extinction,absorption, asymmetry
c           fg1ext=fog1 .2km vis extinction  fg1abs=fog1 absorption
c           fg1sym=fog1 asymmetry factors
c           fg2ext=fog2 .5km vis extinction  fg2abs=fog2 absorption
c           fg2sym=fog2 asymmetry factors
c         >2-10km
c           troext=tropospher extinction  troabs=tropospher absorption
c           trosym=tropospheric asymmetry factors
c         >10-30km
c           bstext=background stratospheric extinction
c           bstabs=background stratospheric absorption
c           bstsym=background stratospheric asymmetry factors
c           avoext=aged volcanic extinction
c           avoabs=aged volcanic absorption
c           avosym=aged volcanic asymmetry factors
c           fvoext=fresh volcanic extinction
c           fvoabs=fresh volcanic absorption
c           fvosym=fresh volcanic asymmetry factors
c         >30-100km
c           dmeext=meteoric dust extinction
c           dmeabs=meteoric dust absorption
c           dmesym=meteoric dust asymmetry factors
c
c     aerosol extinction and absorption data
c
c     modified to include asymmetry data - jan 1986 (a.e.r. inc.)
c
      data rure /
     & 2.09291, 1.74582, 1.60500, 1.00000,  .75203,  .41943,  .24070,
     &  .14709,  .13304,  .12234,  .13247,  .11196,  .10437,  .09956,
     &  .09190,  .08449,  .07861,  .07025,  .07089,  .07196,  .07791,
     &  .04481,  .04399,  .12184,  .12658,  .12829,  .09152,  .08076,
     &  .07456,  .06880,  .06032,  .04949,  .05854,  .06000,  .06962,
     &  .05722,  .06051,  .05177,  .04589,  .04304,
     &  .03582,  .03155,  .02018,  .01469,  .00798,  .00551, 0.,

     & 2.09544, 1.74165, 1.59981, 1.00000,  .75316,  .42171,  .24323,
     &  .15108,  .13608,  .12430,  .13222,  .13823,  .11076,  .10323,
     &  .09475,  .08728,  .08076,  .07639,  .07797,  .07576,  .07943,
     &  .04899,  .04525,  .12165,  .12741,  .12778,  .09032,  .07962,
     &  .07380,  .06880,  .06329,  .05791,  .06646,  .06639,  .07443,
     &  .06304,  .06443,  .05538,  .04867,  .04519,
     &  .03821,  .03374,  .02173,  .01587,  .00862,  .00594, 0.,

     & 2.07082, 1.71456, 1.57962, 1.00000,  .76095,  .43228,  .25348,
     &  .16456,  .14677,  .13234,  .13405,  .20316,  .12873,  .11506,
     &  .10481,  .09709,  .08918,  .09380,  .09709,  .08791,  .08601,
     &  .06247,  .05601,  .11905,  .12595,  .12348,  .08741,  .07703,
     &  .07266,  .07044,  .07443,  .08146,  .08810,  .08563,  .08962,
     &  .08051,  .07677,  .06658,  .05747,  .05184,
     &  .04572,  .04074,  .02689,  .01981,  .01084,  .00714, 0.,

     & 1.66076, 1.47886, 1.40139, 1.00000,  .80652,  .50595,  .32259,
     &  .23468,  .20772,  .18532,  .17348,  .35114,  .20006,  .17386,
     &  .16139,  .15424,  .14557,  .16215,  .16766,  .14994,  .14032,
     &  .12968,  .12601,  .13551,  .13582,  .13228,  .11070,  .09994,
     &  .09873,  .10418,  .13241,  .15924,  .16139,  .15949,  .15778,
     &  .15184,  .13848,  .12563,  .11076,  .09601,
     &  .09312,  .08720,  .06644,  .05264,  .03181,  .02196, 0.0/

      data rura /
     &  .67196,  .11937,  .08506,  .05930,  .05152,  .05816,  .05006,
     &  .01968,  .02070,  .02101,  .05652,  .02785,  .01316,  .00867,
     &  .01462,  .01310,  .01627,  .02013,  .02165,  .02367,  .03538,
     &  .02823,  .03962,  .06778,  .07285,  .08120,  .04032,  .03177,
     &  .02557,  .02342,  .02177,  .02627,  .03943,  .03114,  .03696,
     &  .02956,  .03500,  .03241,  .03297,  .03380,
     &  .03170,  .02794,  .01769,  .01305,  .00730,  .00518, 0.0,

     &  .62968,  .10816,  .07671,  .05380,  .04684,  .05335,  .04614,
     &  .01829,  .01899,  .01962,  .05525,  .06816,  .01652,  .00867,
     &  .01544,  .01373,  .01627,  .02892,  .02829,  .02532,  .03487,
     &  .02835,  .03854,  .06684,  .07272,  .08038,  .03987,  .03247,
     &  .02816,  .02816,  .03101,  .03741,  .04829,  .04032,  .04399,
     &  .03734,  .03956,  .03601,  .03525,  .03563,
     & .03357,  .02965,  .01887,  .01395,  .00782,  .00555, 0.0,

     &  .51899,  .08278,  .05816,  .04082,  .03570,  .04158,  .03620,
     &  .01513,  .01481,  .01633,  .05278,  .13690,  .02494,  .00886,
     &  .01804,  .01582,  .01677,  .04816,  .04367,  .03013,  .03443,
     &  .02930,  .03677,  .06209,  .06911,  .07475,  .03892,  .03494,
     &  .03513,  .03968,  .05152,  .06241,  .06937,  .06203,  .06215,
     &  .05614,  .05209,  .04608,  .04196,  .04095,
     &  .03916,  .03486,  .02262,  .01686,  .00951,  .00674, 0.0,

     &  .21943,  .02848,  .01943,  .01342,  .01171,  .01437,  .01323,
     &  .01152,  .00696,  .01329,  .06108,  .24690,  .05323,  .01430,
     &  .03361,  .02949,  .02652,  .09437,  .08506,  .05348,  .04627,
     &  .04380,  .04557,  .05380,  .05715,  .05899,  .04861,  .05253,
     &  .06171,  .07437,  .10152,  .12019,  .12190,  .11734,  .11411,
     &  .10766,  .09487,  .08430,  .07348,  .06861,
     &  .06936,  .06458,  .04735,  .03761,  .02313,  .01668, 0.0/

      data rurg /
     &  .7581,   .6785,   .6712,   .6479,   .6342,   .6176,   .6334,
     &  .7063,   .7271,   .7463,   .7788,   .7707,   .7424,   .7312,
     &  .7442,   .7516,   .7662,   .7940,   .7886,   .7797,   .7664,
     &  .8525,   .8700,   .5846,   .5570,   .5992,   .6159,   .6271,
     &  .6257,   .6374,   .6546,   .6861,   .6859,   .6120,   .5570,
     &  .5813,   .5341,   .5284,   .5137,   .4348,   .4223,   .3775,
     &  .3435,   .3182,   .2791,   .2494,   .0000,

     &  .7632,   .6928,   .6865,   .6638,   .6498,   .6314,   .6440,
     &  .7098,   .7303,   .7522,   .7903,   .7804,   .7380,   .7319,
     &  .7508,   .7584,   .7738,   .8071,   .7929,   .7843,   .7747,
     &  .8507,   .8750,   .6112,   .5851,   .6272,   .6466,   .6616,
     &  .6653,   .6798,   .6965,   .7026,   .6960,   .6360,   .5848,
     &  .6033,   .5547,   .5445,   .5274,   .4518,   .4318,   .3863,
     &  .3516,   .3257,   .2853,   .2548,   .0000,

     &  .7725,   .7240,   .7197,   .6997,   .6858,   .6650,   .6702,
     &  .7181,   .7378,   .7653,   .8168,   .7661,   .7286,   .7336,
     &  .7654,   .7735,   .7910,   .8303,   .8025,   .7957,   .7946,
     &  .8468,   .8734,   .6831,   .6619,   .6994,   .7250,   .7449,
     &  .7547,   .7665,   .7644,   .7265,   .7170,   .6769,   .6409,
     &  .6442,   .6031,   .5854,   .5646,   .4977,   .4602,   .4127,
     &  .3751,   .3476,   .3048,   .2721,   .0000,

     &  .7778,   .7793,   .7786,   .7717,   .7628,   .7444,   .7365,
     &  .7491,   .7609,   .7921,   .8688,   .7537,   .7294,   .7413,
     &  .7928,   .8016,   .8225,   .8761,   .8359,   .8285,   .8385,
     &  .8559,   .8654,   .8414,   .8415,   .8527,   .8740,   .8903,
     &  .8952,   .8923,   .8611,   .8033,   .7989,   .7758,   .7632,
     &  .7508,   .7314,   .7091,   .6867,   .6419,   .5790,   .5259,
     &  .4749,   .4415,   .3886,   .3489,   .0000/

      data urbe /
     & 1.88816, 1.63316, 1.51867, 1.00000,  .77785,  .47095,  .30006,
     &  .21392,  .19405,  .17886,  .18127,  .16133,  .14785,  .14000,
     &  .12715,  .11880,  .11234,  .10601,  .10500,  .10361,  .10342,
     &  .08766,  .08652,  .11937,  .12139,  .12297,  .09797,  .09057,
     &  .08595,  .08196,  .07563,  .06696,  .07209,  .06842,  .07177,
     &  .06354,  .06177,  .05373,  .04728,  .04051,
     &  .03154,  .02771,  .01759,  .01278,  .00693,  .00480, 0.0,

     & 1.95582, 1.64994, 1.53070, 1.00000,  .77614,  .46639,  .29487,
     &  .21051,  .18943,  .17285,  .17209,  .21418,  .15354,  .14051,
     &  .12728,  .11861,  .11089,  .11329,  .11323,  .10563,  .10247,
     &  .08696,  .08361,  .12013,  .12418,  .12304,  .09614,  .08842,
     &  .08487,  .08285,  .08361,  .08430,  .08880,  .08449,  .08601,
     &  .07835,  .07323,  .06367,  .05500,  .04747,
     &  .03901,  .03454,  .02240,  .01638,  .00891,  .00612, 0.0,

     & 1.96430, 1.64032, 1.52392, 1.00000,  .77709,  .46253,  .28690,
     &  .20310,  .17981,  .16101,  .15614,  .26475,  .15456,  .13563,
     &  .12215,  .11361,  .10500,  .11715,  .11753,  .10392,  .09766,
     &  .08443,  .08057,  .10943,  .11342,  .11063,  .08703,  .08025,
     &  .07886,  .08032,  .09101,  .10070,  .10386,  .09943,  .09886,
     &  .09152,  .08247,  .07152,  .06089,  .05253,
     &  .04582,  .04091,  .02717,  .02008,  .01103,  .00754, 0.0,

     & 1.41266, 1.33816, 1.29114, 1.00000,  .83646,  .55025,  .35342,
     &  .25285,  .21576,  .18310,  .16215,  .37854,  .20494,  .16665,
     &  .14778,  .13892,  .12943,  .15525,  .15709,  .13513,  .12481,
     &  .11759,  .11494,  .11487,  .11329,  .11108,  .09911,  .09209,
     &  .09342,  .10120,  .13177,  .15696,  .15766,  .15513,  .15203,
     &  .14532,  .13038,  .11785,  .10411,  .09101,
     &  .08907,  .08399,  .06579,  .05337,  .03372,  .02379, 0.0/

      data urba /
     &  .78437,  .58975,  .54285,  .36184,  .29222,  .20886,  .15658,
     &  .12329,  .11462,  .10747,  .11797,  .10025,  .08759,  .08184,
     &  .07506,  .07006,  .06741,  .06601,  .06544,  .06449,  .06665,
     &  .06278,  .06949,  .07316,  .07462,  .08101,  .05753,  .05272,
     &  .04899,  .04734,  .04494,  .04443,  .05133,  .04348,  .04443,
     &  .03994,  .03981,  .03633,  .03468,  .03146,
     &  .02809,  .02471,  .01556,  .01145,  .00639,  .00454, 0.0,

     &  .69032,  .49367,  .45165,  .29741,  .24070,  .17399,  .13146,
     &  .10354,  .09589,  .09025,  .10411,  .15101,  .07880,  .06949,
     &  .06570,  .06095,  .05829,  .07171,  .06797,  .05975,  .06013,
     &  .05589,  .06051,  .07139,  .07494,  .07956,  .05525,  .05184,
     &  .05089,  .05291,  .05886,  .06380,  .06880,  .06127,  .06019,
     &  .05525,  .05070,  .04500,  .04076,  .03741,
     &  .03400,  .03010,  .01926,  .01427,  .00800,  .00567, 0.0,

     &  .54848,  .37101,  .33734,  .21949,  .17785,  .12968,  .09854,
     &  .07804,  .07165,  .06791,  .08563,  .19639,  .06722,  .05316,
     &  .05316,  .04886,  .04620,  .07570,  .06899,  .05291,  .05101,
     &  .04734,  .05025,  .06171,  .06570,  .06854,  .04892,  .04797,
     &  .05057,  .05665,  .07127,  .08095,  .08411,  .07728,  .07475,
     &  .06886,  .06019,  .05222,  .04538,  .04171,
     &  .03911,  .03486,  .02271,  .01697,  .00961,  .00681, 0.0,

     &  .15975,  .10000,  .09013,  .05785,  .04671,  .03424,  .02633,
     &  .02525,  .01975,  .02354,  .06241,  .26690,  .05810,  .02285,
     &  .03810,  .03386,  .03044,  .09627,  .08557,  .05405,  .04576,
     &  .04392,  .04424,  .04671,  .04791,  .04861,  .04684,  .05177,
     &  .06158,  .07475,  .10342,  .12146,  .12177,  .11734,  .11335,
     &  .10608,  .09171,  .08063,  .06968,  .06475,
     &  .06559,  .06131,  .04591,  .03714,  .02365,  .01734, 0.0/

      data urbg /
     &  .7785,   .7182,   .7067,   .6617,   .6413,   .6166,   .6287,
     &  .6883,   .7070,   .7243,   .7370,   .7446,   .7391,   .7371,
     &  .7414,   .7435,   .7466,   .7543,   .7498,   .7424,   .7270,
     &  .7674,   .7850,   .5880,   .5616,   .5901,   .6159,   .6238,
     &  .6240,   .6281,   .6306,   .6298,   .6252,   .5785,   .5378,
     &  .5512,   .5072,   .4930,   .4709,   .4009,   .4110,   .3672,
     &  .3344,   .3093,   .2717,   .2426,   .0000,

     &  .7906,   .7476,   .7385,   .6998,   .6803,   .6536,   .6590,
     &  .7066,   .7258,   .7484,   .7769,   .7405,   .7351,   .7459,
     &  .7625,   .7673,   .7759,   .7910,   .7732,   .7703,   .7644,
     &  .7966,   .8142,   .6635,   .6428,   .6700,   .6935,   .7050,
     &  .7092,   .7145,   .7094,   .6762,   .6684,   .6316,   .5997,
     &  .6013,   .5625,   .5433,   .5198,   .4552,   .4387,   .3928,
     &  .3575,   .3310,   .2899,   .2588,   .0000,

     &  .7949,   .7713,   .7650,   .7342,   .7162,   .6873,   .6820,
     &  .7131,   .7312,   .7583,   .8030,   .7171,   .7185,   .7400,
     &  .7698,   .7778,   .7923,   .8142,   .7864,   .7867,   .7891,
     &  .8147,   .8298,   .7276,   .7136,   .7361,   .7590,   .7729,
     &  .7783,   .7808,   .7624,   .7094,   .7022,   .6714,   .6480,
     &  .6417,   .6104,   .5887,   .5651,   .5058,   .4692,   .4212,
     &  .3825,   .3549,   .3112,   .2778,   .0000,

     &  .7814,   .7993,   .7995,   .7948,   .7870,   .7682,   .7751,
     &  .7501,   .7565,   .7809,   .8516,   .7137,   .7039,   .7241,
     &  .7728,   .7846,   .8093,   .8576,   .8125,   .8140,   .8304,
     &  .8472,   .8549,   .8525,   .8569,   .8640,   .8853,   .9017,
     &  .9061,   .9021,   .8685,   .8126,   .8091,   .7897,   .7802,
     &  .7691,   .7550,   .7353,   .7146,   .6754,   .6134,   .5601,
     &  .5056,   .4701,   .4134,   .3714,   .0000/

      data ocne /
     & 1.47576, 1.32614, 1.26171, 1.00000,  .88133,  .70297,  .56487,
     &  .46006,  .42044,  .38310,  .35076,  .42266,  .32278,  .28810,
     &  .24905,  .21184,  .16734,  .14791,  .21532,  .15076,  .12057,
     &  .10038,  .10703,  .15070,  .15665,  .14639,  .10228,  .08367,
     &  .07373,  .06829,  .05044,  .04373,  .04962,  .06158,  .07703,
     &  .07234,  .06297,  .05481,  .05329,  .08741,
     &  .04608,  .03959,  .02382,  .01712,  .00936,  .00665, 0.0,

     & 1.36924, 1.25443, 1.20835, 1.00000,  .91367,  .77089,  .64987,
     &  .54886,  .50247,  .45038,  .38209,  .50589,  .43766,  .38076,
     &  .31658,  .27475,  .22215,  .21019,  .27570,  .21057,  .16949,
     &  .14209,  .14215,  .16956,  .17082,  .16025,  .11665,  .09759,
     &  .09215,  .09373,  .10532,  .12570,  .13000,  .13633,  .14291,
     &  .13506,  .11475,  .09658,  .08291,  .10348,
     &  .06693,  .05786,  .03522,  .02519,  .01358,  .00954, 0.0,

     & 1.22259, 1.14627, 1.11842, 1.00000,  .94766,  .87538,  .80418,
     &  .72930,  .68582,  .62165,  .49962,  .67949,  .66468,  .59253,
     &  .49551,  .44671,  .37886,  .35924,  .43367,  .37019,  .30842,
     &  .26437,  .25228,  .24905,  .23975,  .22766,  .17804,  .15316,
     &  .15373,  .16791,  .22361,  .28348,  .28677,  .29082,  .29038,
     &  .27810,  .23867,  .20209,  .16430,  .14943,
     &  .12693,  .11177,  .07095,  .05084,  .02690,  .01838, 0.0,

     & 1.09133, 1.06601, 1.05620, 1.00000,  .97506,  .94791,  .94203,
     &  .93671,  .92867,  .90411,  .80253,  .89222,  .94462,  .92146,
     &  .85797,  .82595,  .76747,  .68646,  .78209,  .75266,  .68658,
     &  .62722,  .60228,  .56335,  .53728,  .51861,  .43449,  .37196,
     &  .35899,  .37316,  .46854,  .58234,  .58690,  .60348,  .60563,
     &  .60000,  .55392,  .50367,  .43576,  .35949,
     &  .34729,  .32254,  .23600,  .17953,  .10071,  .06714, 0.0/

      data ocna /
     &  .30987,  .04354,  .02880,  .01797,  .01468,  .01766,  .01582,
     &  .00816,  .01146,  .01677,  .03310,  .03380,  .00715,  .00443,
     &  .00500,  .00601,  .00753,  .01595,  .02943,  .00994,  .01367,
     &  .01671,  .02538,  .03481,  .03405,  .03601,  .01608,  .01310,
     &  .01152,  .01082,  .01070,  .01563,  .02063,  .03171,  .03810,
     &  .03741,  .03804,  .03759,  .04209,  .07892,
     &  .04347,  .03754,  .02269,  .01649,  .00917,  .00657, 0.0,

     &  .23367,  .03127,  .02070,  .01297,  .01063,  .01285,  .01190,
     &  .00937,  .00911,  .01576,  .05576,  .23487,  .03949,  .00905,
     &  .02057,  .01816,  .01665,  .08025,  .08044,  .03677,  .03139,
     &  .03190,  .03766,  .04532,  .04544,  .04715,  .03405,  .03614,
     &  .04329,  .05424,  .07823,  .09728,  .10057,  .10247,  .10222,
     &  .09551,  .08241,  .07158,  .06506,  .09203,
     &  .06133,  .05332,  .03258,  .02366,  .01308,  .00932, 0.0,

     &  .13025,  .01557,  .01013,  .00646,  .00532,  .00665,  .00722,
     &  .01335,  .00728,  .01810,  .09835,  .37329,  .09703,  .01968,
     &  .05114,  .04342,  .03709,  .17456,  .16468,  .08785,  .06880,
     &  .06589,  .06791,  .07247,  .07329,  .07449,  .07025,  .07962,
     &  .09899,  .12481,  .17867,  .22019,  .22228,  .22051,  .21595,
     &  .20335,  .17278,  .14677,  .12171,  .12430,
     &  .10890,  .09644,  .06106,  .04465,  .02457,  .01732, 0.0,

     &  .03506,  .00323,  .00215,  .00139,  .00114,  .00171,  .00532,
     &  .03082,  .01101,  .03741,  .20101,  .47608,  .21165,  .05234,
     &  .12886,  .11215,  .09684,  .32810,  .31778,  .20513,  .16658,
     &  .15956,  .15842,  .15905,  .15968,  .16051,  .16506,  .18323,
     &  .21709,  .25652,  .33222,  .39639,  .39854,  .40297,  .40025,
     &  .39025,  .35468,  .32006,  .27715,  .25348,
     &  .25632,  .23876,  .17092,  .13198,  .07692,  .05407, 0.0/

      data ocng /
     &  .7516,   .6960,   .6920,   .6756,   .6767,   .6844,   .6936,
     &  .7055,   .7110,   .7177,   .7367,   .6287,   .6779,   .6784,
     &  .6599,   .6659,   .6859,   .6887,   .6095,   .6558,   .6665,
     &  .6697,   .6594,   .5851,   .5644,   .5760,   .5903,   .5991,
     &  .6024,   .5979,   .6087,   .5837,   .5763,   .5348,   .4955,
     &  .4821,   .4635,   .4373,   .3944,   .2344,   .2754,   .2447,
     &  .2266,   .2088,   .1766,   .1481,   .0000,

     &  .7708,   .7288,   .7243,   .7214,   .7211,   .7330,   .7445,
     &  .7579,   .7649,   .7790,   .8182,   .7673,   .7171,   .7205,
     &  .7235,   .7251,   .7397,   .7537,   .6934,   .7137,   .7193,
     &  .7206,   .7151,   .6732,   .6620,   .6696,   .6821,   .6895,
     &  .6898,   .6819,   .6556,   .5925,   .5869,   .5511,   .5284,
     &  .5124,   .4912,   .4646,   .4302,   .3124,   .3101,   .2752,
     &  .2529,   .2335,   .2021,   .1738,   .0000,

     &  .7954,   .7782,   .7752,   .7717,   .7721,   .7777,   .7872,
     &  .8013,   .8089,   .8301,   .8844,   .8332,   .7557,   .7597,
     &  .7823,   .7822,   .7944,   .8157,   .7712,   .7738,   .7784,
     &  .7807,   .7800,   .7682,   .7659,   .7692,   .7780,   .7828,
     &  .7776,   .7621,   .7115,   .6342,   .6294,   .5999,   .5854,
     &  .5700,   .5512,   .5265,   .4996,   .4236,   .3765,   .3357,
     &  .3066,   .2830,   .2466,   .2184,   .0000,

     &  .8208,   .8270,   .8260,   .8196,   .8176,   .8096,   .8096,
     &  .8202,   .8255,   .8520,   .9228,   .8950,   .7965,   .7847,
     &  .8242,   .8244,   .8376,   .8857,   .8463,   .8332,   .8379,
     &  .8441,   .8467,   .8502,   .8534,   .8562,   .8688,   .8789,
     &  .8785,   .8683,   .8252,   .7562,   .7519,   .7261,   .7141,
     &  .6980,   .6789,   .6540,   .6294,   .5783,   .5100,   .4595,
     &  .4164,   .3868,   .3404,   .3042,   .0000/

      data troe /
     & 2.21222, 1.82753, 1.67032, 1.00000,  .72424,  .35272,  .15234,
     &  .05165,  .03861,  .02994,  .04671,  .02462,  .01538,  .01146,
     &  .01032,  .00816,  .00861,  .00994,  .01057,  .01139,  .01747,
     &  .01494,  .02418,  .03165,  .03386,  .04247,  .01601,  .01215,
     &  .00937,  .00861,  .00823,  .01139,  .01924,  .01234,  .01348,
     &  .01114,  .01297,  .01266,  .01418,  .01487,
     &  .01543,  .01321,  .00793,  .00582,  .00330,  .00239, 0.0,

     & 2.21519, 1.82266, 1.66557, 1.00000,  .72525,  .35481,  .15449,
     &  .05475,  .04044,  .03082,  .04620,  .05272,  .01867,  .01266,
     &  .01127,  .00886,  .00886,  .01449,  .01399,  .01228,  .01728,
     &  .01475,  .02285,  .03215,  .03494,  .04285,  .01652,  .01304,
     &  .01101,  .01120,  .01297,  .01753,  .02468,  .01741,  .01766,
     &  .01513,  .01557,  .01456,  .01532,  .01582,
     &  .01619,  .01386,  .00832,  .00610,  .00346,  .00251, 0.0,

     & 2.19082, 1.79462, 1.64456, 1.00000,  .73297,  .36443,  .16278,
     &  .06468,  .04658,  .03399,  .04538,  .11892,  .02835,  .01646,
     &  .01386,  .01076,  .00968,  .02551,  .02222,  .01468,  .01690,
     &  .01437,  .01994,  .03127,  .03513,  .04076,  .01722,  .01513,
     &  .01519,  .01791,  .02538,  .03272,  .03816,  .03038,  .02886,
     &  .02551,  .02228,  .01937,  .01804,  .01791,
     &  .01798,  .01539,  .00924,  .00678,  .00384,  .00278, 0.0,

     & 1.75696, 1.54829, 1.45962, 1.00000,  .77816,  .43139,  .21778,
     &  .11329,  .08101,  .05506,  .04943,  .25291,  .06816,  .03703,
     &  .02601,  .01968,  .01468,  .04962,  .04247,  .02234,  .01797,
     &  .01532,  .01633,  .02259,  .02487,  .02595,  .01728,  .01892,
     &  .02399,  .03247,  .05285,  .06462,  .06608,  .05930,  .05525,
     &  .04861,  .03753,  .02968,  .02348,  .02165,
     &  .02152,  .01841,  .01104,  .00809,  .00458,  .00332, 0.0/

      data troa /
     &  .69671,  .09905,  .06563,  .04101,  .03354,  .03627,  .02810,
     &  .00873,  .00918,  .00930,  .03215,  .01285,  .00513,  .00316,
     &  .00557,  .00494,  .00646,  .00867,  .00937,  .01025,  .01646,
     &  .01481,  .02418,  .02886,  .03070,  .04032,  .01494,  .01139,
     &  .00873,  .00816,  .00797,  .01133,  .01911,  .01215,  .01329,
     &  .01101,  .01291,  .01266,  .01418,  .01487,
     &  .01543,  .01321,  .00793,  .00582,  .00330,  .00239, 0.0,

     &  .65000,  .08791,  .05816,  .03652,  .02994,  .03278,  .02557,
     &  .00810,  .00842,  .00867,  .03139,  .03949,  .00646,  .00316,
     &  .00595,  .00519,  .00646,  .01304,  .01247,  .01095,  .01620,
     &  .01449,  .02278,  .02930,  .03184,  .04063,  .01544,  .01234,
     &  .01044,  .01076,  .01272,  .01741,  .02462,  .01722,  .01747,
     &  .01506,  .01551,  .01456,  .01532,  .01582,
     &  .01619,  .01386,  .00832,  .00610,  .00346,  .00251, 0.0,

     &  .52804,  .06367,  .04158,  .02633,  .02184,  .02443,  .01937,
     &  .00658,  .00646,  .00709,  .02949,  .10013,  .00968,  .00310,
     &  .00677,  .00582,  .00646,  .02361,  .01994,  .01266,  .01544,
     &  .01386,  .01968,  .02848,  .03203,  .03854,  .01620,  .01449,
     &  .01462,  .01747,  .02513,  .03253,  .03797,  .03019,  .02861,
     &  .02538,  .02215,  .01930,  .01797,  .01791,
     &  .01797,  .01539,  .00924,  .00677,  .00384,  .00278, 0.0,

     &  .19829,  .01842,  .01215,  .00791,  .00665,  .00778,  .00652,
     &  .00361,  .00253,  .00399,  .02570,  .20690,  .01715,  .00316,
     &  .00873,  .00728,  .00658,  .04481,  .03525,  .01646,  .01405,
     &  .01310,  .01468,  .01956,  .02184,  .02367,  .01608,  .01816,
     &  .02342,  .03203,  .05234,  .06399,  .06538,  .05867,  .05456,
     &  .04810,  .03715,  .02949,  .02335,  .02158,
     &  .02149,  .01840,  .01104,  .00809,  .00458,  .00332, 0.0/

      data trog /
     &  .7518,   .6710,   .6638,   .6345,   .6152,   .5736,   .5280,
     &  .4949,   .4700,   .4467,   .4204,   .4028,   .3777,   .3563,
     &  .3150,   .2919,   .2695,   .2465,   .2402,   .2313,   .2101,
     &  .1760,   .1532,   .2091,   .2079,   .1843,   .1811,   .1687,
     &  .1626,   .1526,   .1356,   .1030,   .0962,   .1024,   .1086,
     &  .0928,   .0836,   .0643,   .0451,   .0290,   .0156,   .0118,
     &  .0076,   .0050,   .0024,   .0015,   .0000,

     &  .7571,   .6858,   .6790,   .6510,   .6315,   .5887,   .5418,
     &  .5075,   .4829,   .4598,   .4338,   .4043,   .3890,   .3680,
     &  .3259,   .3026,   .2800,   .2541,   .2494,   .2414,   .2196,
     &  .1873,   .1657,   .2123,   .2110,   .1890,   .1836,   .1709,
     &  .1640,   .1534,   .1354,   .1044,   .0984,   .1026,   .1073,
     &  .0935,   .0842,   .0661,   .0477,   .0309,   .0171,   .0129,
     &  .0084,   .0056,   .0027,   .0017,   .0000,

     &  .7667,   .7176,   .7128,   .6879,   .6690,   .6255,   .5769,
     &  .5403,   .5167,   .4947,   .4703,   .4143,   .4190,   .3993,
     &  .3563,   .3325,   .3095,   .2767,   .2751,   .2693,   .2464,
     &  .2175,   .1992,   .2247,   .2215,   .2042,   .1952,   .1814,
     &  .1726,   .1604,   .1398,   .1111,   .1065,   .1068,   .1086,
     &  .0984,   .0888,   .0724,   .0549,   .0358,   .0216,   .0166,
     &  .0109,   .0073,   .0036,   .0023,   .0000,

     &  .7696,   .7719,   .7710,   .7606,   .7478,   .7142,   .6727,
     &  .6381,   .6201,   .6050,   .5912,   .4849,   .5137,   .5019,
     &  .4625,   .4389,   .4169,   .3696,   .3707,   .3708,   .3473,
     &  .3232,   .3112,   .3022,   .2938,   .2850,   .2675,   .2494,
     &  .2347,   .2165,   .1857,   .1536,   .1509,   .1441,   .1416,
     &  .1354,   .1245,   .1088,   .0905,   .0614,   .0440,   .0354,
     &  .0257,   .0179,   .0089,   .0059,   .0000/

c
      if(iaer.le.0.) then
        wlbaer(1)=0
        return
      endif
    
      abaer=0.                       ! angstrom exponent      

      rhum=max(zero,min(one,humid))
      if(rhum.lt.rhzone(2)) then
        j=1
      elseif(rhum.lt.rhzone(3)) then
        j=2
      else
        j=3
      endif
      
      wt=(rhum-rhzone(j))/(rhzone(j+1)-rhzone(j))      
      
      nwlbaer=naerw

      do 20 i=1,naerw
        wlbaer(i)=awl(i)
        if(abs(iaer).eq.1) then
          ex1=max(rure(i,j),tiny)
          ex2=max(rure(i,j+1),tiny)
          ab1=max(rura(i,j),tiny)
          ab2=max(rura(i,j+1),tiny)
          g1=max(rurg(i,j),tiny)
          g2=max(rurg(i,j+1),tiny)
        elseif(abs(iaer).eq.2) then
          ex1=max(urbe(i,j),tiny)
          ex2=max(urbe(i,j+1),tiny)
          ab1=max(urba(i,j),tiny)
          ab2=max(urba(i,j+1),tiny)
          g1=max(urbg(i,j),tiny)
          g2=max(urbg(i,j+1),tiny)
        elseif(abs(iaer).eq.3) then
          ex1=max(ocne(i,j),tiny)
          ex2=max(ocne(i,j+1),tiny)
          ab1=max(ocna(i,j),tiny)
          ab2=max(ocna(i,j+1),tiny)
          g1=max(ocng(i,j),tiny)
          g2=max(ocng(i,j+1),tiny)
        elseif(abs(iaer).eq.4) then
          ex1=max(troe(i,j),tiny)
          ex2=max(troe(i,j+1),tiny)
          ab1=max(troa(i,j),tiny)
          ab2=max(troa(i,j+1),tiny)
          g1=max(trog(i,j),tiny)
          g2=max(trog(i,j+1),tiny)
        endif
        if(abs(iaer).ge.1.and.abs(iaer).le.4) then
          aerext(i)=ex1*(ex2/ex1)**wt
          aerabs(i)=ab1*(ab2/ab1)**wt
          aerasm(i)=g1*(g2/g1)**wt
        endif
 20   continue

c        write(*,1000) 'wlbaer:',wlbaer
c        write(*,1000) 'aerext:',aerext
c        write(*,1000) 'aerabs:',aerabs
c        write(*,1000) 'aerasm:',aerasm
 1000 format(a/5(10es11.3,/))
      return
      end subroutine stdaer
c
c-----------------------------------------------------------------------

      end module aeroblk

c=======================================================================
      function aeroden(zz)

c purpose:   find number density of boundary layer aerosols, aeroden,
c            at a given altitude, zz, and for a specified visibility
c input:
c   zz       altitude   (km)
c
c output:
c            aerosol density at altitude zz
c

      use params, only: zero,one,kr
      use aeroblk, only: dbaer, zbaer, nzbaer
      implicit none

      integer :: k,kp

      real(kr) :: zz,z,f,aeroden
      real(kr), parameter :: z1=0., z2=100.

      z=max(z1,min(z2,zz))
      aeroden=0.
c
      if(z.gt.zbaer(nzbaer)) return
      call locate(zbaer,nzbaer,z,k)
      kp=k+1
      f=(z-zbaer(k))/(zbaer(kp)-zbaer(k))
      
      if(min(dbaer(k),dbaer(kp)).le.zero) then
        aeroden=max(dbaer(k)*(1.-f)+dbaer(kp)*f,0.0_kr)
      else
        aeroden=dbaer(k)*(dbaer(kp)/dbaer(k))**f
      endif
      
        !write(*,*) 'aeroden k,kp,z,dbaer(k),dbaer(kp),f,aeroden'
        !write(*,'(2i5,5es11.3)') k,kp,z,dbaer(k),dbaer(kp),f,aeroden
        
      RETURN
      end function aeroden

c=======================================================================
      subroutine tauaero(wl,nz,z,nosct,dtaua,waer,nmom,pmom,idb)
c
c purpose:  
c     compute the optical depth and scattering parameters at all
c     atmospheric levels of both boundary layer and stratospheric
c     aerosols
c
c input:
c   wl      wavelength
c   nz      number of atmospheric levels
c   z       altitude array (km)
c   idb     if 1 then print diagnostic
c   nosct   0 => normal scattering; 1 => no scattering assumption
C
c output:
c   dtaua   aerosol optical depth increment at each atmospheric layer
c   waer    single scatter albedo due to aerosols at each atmospheric layer
c
c input/output:
c   nmom   number of phase function moments
c   pmom   input:  (moments of phase function of cloud)*taucld*wcld
c          output: (moments of phase function clouds)*taucld*wcld +
c                  (moments of phase function aerosol)*tauaer*waer
c
c NOTE:     visibility is defined as 3.912/extinction_coefficient(km-1)
c           for a horizontal path 
c           (exp(-3.912) = 0.02 = visible contrast threshold)
c
c revisions:
c
c  971028: corrected expression for net asymmetry factor including both 
c          stratospheric and boundary layer aerosols (suggested by Xiang Li)
c
c  040428: corrected bug that prevented computation of stratospheric aerosols
c          when iaer=0
c
      use params, only: mxly, maxmom, zero, kr, zip
      use aeroblk, only: naerz, iaer, imoma, tbaer, jaer, zaer,
     &     taerst, npmaer, aestrat, phaerw, aerbwi

      implicit none

      real(kr), intent(in)    :: wl,z(*)
      real(kr), intent(out)   :: dtaua(*),waer(*)
      real(kr), intent(inout) :: pmom(0:maxmom,*)
      real(kr), save          :: dtsv(mxly)

      real(kr) :: dtauab(mxly),gaer(mxly),pm(0:maxmom),extinc,wa,ga,dt

      integer, intent(in) :: nz,nosct,nmom,idb
      integer             :: i,j,nl,ja,namom 
      integer, external   :: numset
      integer, save       :: init=1, laer(naerz)

      logical, parameter  :: debug=.false.

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                                ! find stratospheric aerosol layers

      if(any(jaer(1:naerz).ne.0).and.init.eq.1) then
        call zlayer(nz,z,numset(zero,taerst,naerz),zaer,laer)
        if(debug) then
          write(*,'(a,5f10.4)') '      zaer:' ,zaer
          write(*,'(a,5i10)')   '      laer:' ,laer
        endif
      endif

      if(iaer.eq.0) then        ! no BL aerosols

        dtaua(1:nz)=0.
        dtauab(1:nz)=0.
        waer(1:nz)=0.
        gaer(1:nz)=0.

      elseif(iaer.eq.-1) then   ! read in BLA model from aerosol.dat

        gaer(1:nz)=pmom(1,1:nz)            ! save cloud effect
        call aeread(imoma,nz,wl,dtauab,waer,pmom)
        gaer(1:nz)=pmom(1,1:nz)-gaer(1:nz) ! subtract cloud effect 

      else                      ! use one of the standard BLA models

c integrate vertical profile of boundary layer aerosol density
c dtsv = sigma*n*dz, where sigma is geometric cross-section
c dtsv is independent of wavelenth and is is saved for later calls to tauaero
      
        if(init.eq.1) call denprfl(nz,z,dtsv)
        
c compute scattering parameters (extinction efficiency , single 
c scattering albedo and asymmetry factor) for boundary layer 
c aerosols at a given wavelength. Note that extinction efficiency
c is normalized to 1 at 550nm.
        
        if (debug) write(*,*) ' call aerbwi'
        call aerbwi(wl,extinc,wa,ga)
        
c nosct options:
c  1.  tauaer is absorption optical depth
c  2.  tauaer is extinction optical depth, ssa = 0
c  3.  tauaer is absorption + backscatter optical depth, ssa = 0

        if(nosct.eq.1) extinc = extinc*(1.-wa)
        if(nosct.eq.3) extinc = extinc*(1.-wa*ga)
        if (nosct.ne.0) then
          wa=0.
          ga=0.
        endif 
        
c use standard phasefunction if imoma is set

        if(imoma.gt.0) then
          namom=nmom
          call getmom(imoma,ga,nmom,pm(0))
        else
          namom=npmaer
          call phaerw(wl,pm(1))
        endif
c        print *,"z,dtauab,extinc,dtsv,wa,pm"
        do i=1,nz
          dtauab(i)=extinc*dtsv(i)
c          print '(6es11.3)', z(nz-i+1),dtauab(i),extinc,dtsv(i),wa,pm(1)
          waer(i)=wa
          gaer(i)=pm(1)*waer(i)*dtauab(i)
          do j=1,namom
            pmom(j,i)=pmom(j,i)+pm(j)*dtauab(i)*waer(i)
          enddo
        enddo
        
      endif                     ! end use one of the standard models

c      print *,'total optical depth ',sum(dtauab(1:nz))

c add in stratospheric aerosols at specified scattering layers

c      print *, "boundary layer"
      do i=1,nz
        dtaua(i)=dtauab(i)
      enddo

      do i=1,naerz

        if(jaer(i).ne.0 .and. taerst(i).gt.0.) then
          nl=laer(i) 
          ja=jaer(i)

c find scattering parameters for stratospheric
c aerosols at a given wavelength

          call aestrat(ja,wl,extinc,wa,ga)
          dt=taerst(i)*extinc
          
          call getmom(3,ga,nmom,pm(0))

          do j=1,nmom
            pmom(j,nl)=pmom(j,nl)+pm(j)*dt*wa
          enddo

          waer(nl)=(waer(nl)*dtaua(nl)+wa*dt)/(dtaua(nl)+dt)
          dtaua(nl)=dtaua(nl)+dt
          gaer(nl)=gaer(nl)+pm(1)*dt*waer(nl)
          
        endif
      enddo

      if (idb.gt.0) then
        where(dtaua(1:nz)*waer(1:nz).ne.0)
     &       gaer(1:nz)=gaer(1:nz)/(dtaua(1:nz)*waer(1:nz))
 
        print '(a3,es11.3,a8,i3,a8,5i4)', 
     &       'wl=', wl,'   iaer=',iaer,'   jaer=',jaer 
        print '(a3,5a11)',
     &       'i','z','dtaua','dtauab','waer','gaer'
        print '(i3,5es11.3)',(i,z(nz-i+1),dtaua(i),dtauab(i),waer(i),
     &       gaer(i),i=1,nz)
        print '(3x,11x,2a11)','---------','----------'
        print '(3x,11x,2es11.3)',sum(dtaua(1:nz)),sum(dtauab(1:nz))
        print *
      end if

      init=0                    ! first call is done

      return

      end
c========================================================================
      subroutine denprfl(nz,z,dtsv)

c  input:
c    nz     number of altitude layers
c    z      altitude array (km)
c
c  output:
c    dtsv  the array of values sigma * n * dz where sigma is the geometric
c          cross-section, n is the aerosol number density.  The optical
c          depth increments of boundary layer aerosols is ext*dtsv where
c          ext is the wavelength dependent extinction efficiency.

      use params, only: mxly, zero, wl55, kr, zip
      use aeroblk, only: iaer, nzbaer, zbaer, ndbaer,
     &     dbaer, vis, tbaer, rhaer, stdaer, usraer, aerzstd, aerbwi

      implicit none

      real(kr), intent(in)    :: z(*)
      real(kr), intent(out)   :: dtsv(*)
      real(kr), parameter     :: visfac=3.912
      real(kr), external      :: aeroden
      real(kr) :: sigma

      real(kr) :: ext55,w55,g55,vint,aden0

      integer, intent(in) :: nz
      integer             :: i
      integer, external   :: numset

      logical, parameter  :: debug=.false.
      
      sigma=0.
      
      nzbaer=numset(zip,zbaer,mxly)
      ndbaer=numset(zip,dbaer,mxly)
      
      if(ndbaer.gt.0) then
        if(nzbaer.eq.1) then
          write(*,*) 'Error -- only one value of zbaer set'
          stop
        endif
        if(nzbaer.ne.ndbaer.and.nzbaer.gt.1) then
          write(*,*) 'Error -- number of elements must match:'
          write(*,'(a,/,(10es11.3))') 'zbaer',(zbaer(i),i=1,nzbaer)
          write(*,'(a,/,(10es11.3))') 'dbaer',(dbaer(i),i=1,ndbaer)
          stop
        endif
        if(nzbaer.eq.0) then
          nzbaer=ndbaer
          zbaer=z(1:nzbaer)
        endif
      else
        call aerzstd            ! copy standard profile to zbaer,dbaer
      endif
      
      if(iaer.eq.5) then
        call usraer(ext55)      ! user defined BLA spectral profile
        if(vis.eq.zip.and.tbaer.eq.zip) tbaer=ext55
      else
        call stdaer(rhaer)      ! standard BLA spectral profile
        call aerbwi(wl55,ext55,w55,g55)
        if(vis.eq.zip.and.tbaer.eq.zip) then
          print *,'must specify either tbaer or vis'
          stop
        endif
      endif
      
      if (debug) write(*,*) ' call aervint',vis,nz
      call aervint(nz,z,dtsv)  ! at this point dtsv = n(i) dz 
      if (debug) write(*,*) ' after call aervint'
      
      if (debug) write(*,*) ' call aerbwi'
      if(ext55.gt.0.) then 
        if(tbaer.ge.0) then
          vint=sum(dtsv(1:nz))
          if(vint.ne.0) sigma=tbaer/(ext55*vint)
        else
          if (debug) write(*,*) ' call aeroden'
          aden0=aeroden(zero)
          sigma=visfac/(ext55*vis*aden0) 
        endif
      endif
      dtsv(1:nz)=sigma*dtsv(1:nz)
      
      end subroutine denprfl
c=======================================================================
 
      subroutine aervint(nz,z,vint)
c
c purpose  
c   compute aerosol density vertical integral (actually the 
c   integral increment at each level, vint=n(i)*dz  (km-1)
c   (the aerosol density can be considered dimensionless since
c   any constant factors are normalized away in tauaero)
c
c input
c   z        monitonically increasing or decreasing array of altitudes (km)
c   nz       number of altitudes
c
c output
c   vint   = n(i) dz, where n(i) is the aerosol density at level i
c            averaged over the interval z(i) to z(i+1)
c            i.e., the column depth increment of aerosol density


      use params, only: kr
      implicit none
      real(kr), parameter :: ztop=100.
      real(kr) :: z(*),vint(*),zu,zd,aden2,aden,dz
      real(kr), external :: aeroden
      integer :: ii,nz,i
      logical :: bup

      bup=z(1).lt.z(nz) ! true if z(1) is at surface
      !print *,'z1,znz,bup ', z1,z(nz),bup

      aden2=aeroden(ztop)
      vint(1)=aden2*5.  ! (aerosol density at nz) * (5 km scale height)

      !print '(6a12)','z(i)','dz','aden','zu','zd','vint(i)'
      !print '(6f12.4)',z(nz),0,aden2,0.,0.,vint(1)

      zu=max(z(1),z(nz))
      do 30 i=2,nz
        ii=i
        if(bup) ii=nz-i+1
        zd=z(ii)
        aden=aeroden(zd)
        dz=zu-zd
        vint(i)=dz*aden
        !print '(6f12.4)',z(ii),dz,aden,zu,zd,vint(i)
        zu=zd
 30   CONTINUE

      return
      end
c===========================================================================
      function relhum(t,h2o)
c
c purpose:
c     return relative humidity for a give temperature and water vapor 
c     density
c
c input:
c   t          temperature (kelvin)
c   h2o        water vapor density (g/m3)
c output:
c   relhum     relative humidity 
c
      use params, only: kr, tzero
      implicit none
      real(kr) :: h2osat, a, t, h2o, relhum

      a=tzero/t
c
c h2osat is the mass density (g/m3) of water vapor at 100% saturation 
c (source: handbook of chemistry and physics, h2o vapor pressure table d-112,
c assuming density related to pressure by ideal gas law)
c
      h2osat=a*exp(18.916758_kr-a*(14.845878_kr+a*2.4918766_kr))
      relhum=h2o/h2osat
      return
      end
c======================================================================
      subroutine aeread(imoma,nz,wl,dtau,wbaer,pmom)
c
c  purpose:  open and read aerosol.dat, perform wavelength interpolation
c
c  input:    
c    imoma   phase function selector
c    nz      number of atmospheric levels 
c    wl      wavelength at which aerosol information is required
c    pmom    (legendre moments, cloud)*tcld*wcld
c
c  output:
c    dtau    array of  aerosol optical depth increments
c    wbaer   array of aerosol single scattering albedo
c    pmom    (legendre moments, cloud)*tcld*wcld +
c            (legendre moments, aerosols)*tauaer*waer
c            
c  method:
c    input file aerosol.dat is opened and read each time wl extends
c    beyond the range spanned by wlbaer(ind) - wlbaer(3-ind).  ind
c    toggles between values 1 and 2.  If aerosol.dat contains
c    information for only one wavelength then that info is assumed
c    spectrally uniform and will be used for all wavelengths.  If
c    aerosol.dat contains info for more than one wavelength,
c    interpolations are done if wl is within the wavelength
c    range. Otherwise if wl is outside of the range, values at the range
c    limits are used for extrapolation. It is best to call aeread with 
c    wavelength arguments that increase with each new call of aeread.
c    However, if aeread is called with a wavelength less than the
c    wavelength range established by previous calls it will rewind
c    aerosol.dat and re-establish the correct wavelengths in
c    wlbaer(ind) and wlbaer(3-ind).
c 
c  file format:
c
c    nn  nmom
c    wl_1
c    (taer_1(i),waer_1(i),(pm_1(j,i),j=1,nmom),i=ns,nz)          
c    wl_2
c    (taer_2(i),waer_2(i),(pm_2(j,i),j=1,nmom),i=ns,nz)
c    .
c    .
c    wl_n
c    (taer_n(i),waer_n(i),(pm_n(j,i),j=1,nmom),i=ns,nz)
c    ...
c    etc.....
c
c    where, wl is the wavelength, with wl_1 < wl_2 < wl_n
c           taer_w (i) is the total optical depth at wavelength w. level i
c           waer_w is the single scattering albedo
c           pm_w(j,i) is the jth moment of the phase function
c           ns=nz-nn+1, 
c           nz is the total number of computational layers
c           nn is the number of aerosol layers read from aerosol.dat
c           
c           taer(nz) is aerosol optical depth at surface
c    

      use params, only: mxly, zero, one, kr, maxmom

      implicit none
      
      real(kr), parameter :: taumin=1.e-30

      real(kr) ::
     &     taer(mxly,2)=0.,
     &     waer(mxly,2)=0.,
     &     gaer(maxmom,mxly,2)=0.,
     &     wlbaer(2)=0.,
     &     wl0=0.

      real(kr) ::
     &     dtau(*), wbaer(*), pmom(0:maxmom,*),
     &     pm(0:maxmom), wl, wt, wlmx, wlinp, gg

      integer :: ind=1,more=1

      integer :: i,nz,nn,ns,nmom,nmomg,imoma

      logical, parameter :: debug=.false.


c mfl commented out
c     data wlbaer/0.,0./
      save taer,waer,gaer,ind,more,wlbaer,ns,nmom,wl0

      if(wlbaer(1).eq.0.) then 
        open(unit=15,file='aerosol.dat',form='formatted',status='old')

        read(15,*) nn, nmom
        ns=nz-nn+1
        if(ns.le.0) then
          print *,'nz  nn ',nz,nn
          print *,'too many layers specified in aerosol.dat'
          stop
        endif

        read(15,*,end=30) wlbaer(1)
        wl0=wlbaer(1)
        do i=ns,nz
          read(15,*,end=20) taer(i,1),waer(i,1),gaer(1:nmom,i,1)
        enddo
        ind=1
      elseif(wl.lt.minval(wlbaer).and.wl.gt.wl0) then
        ! reposition to earlier record if wl < current wavelength interval
        rewind 15
        read(15,*)
        read(15,*) wlbaer(1)
        do i=ns,nz
          taer(i,1)=0.
          waer(i,1)=0.
          gaer(1:nmom,i,1)=0.
          read(15,*,end=20) taer(i,1),waer(i,1),gaer(1:nmom,i,1)
        enddo
        wlbaer(2)=0.
        ind=1
        more=1
      endif
  
      if(more.eq.1) then
        
        do while (wl.gt.maxval(wlbaer))
          more=0
          read(15,*,end=10) wlinp
          ind=3-ind
          do i=ns,nz
            taer(i,ind)=0.
            waer(i,ind)=0.
            gaer(1:nmom,i,ind)=0.
            read(15,*,end=20) taer(i,ind),waer(i,ind),gaer(1:nmom,i,ind)
          enddo
          wlbaer(ind)=wlinp
          if(debug) print '(a,i5,4f12.5)','ind,wl,wbaer ',
     &         ind,wl,wlbaer(1:2)
        enddo
        more=1                  ! didnt hit eof, so there must be more
 10     continue
      endif

      if(wlbaer(2).eq.0) then  
        ! only one wavelength in file, assume spectrally uniform
        wlbaer(1)=.5*wl
        wlbaer(2)=2*wl
        taer(ns:nz,2)=taer(ns:nz,1)
        waer(ns:nz,2)=waer(ns:nz,1)
        gaer(1:nmom,ns:nz,2)=gaer(1:nmom,ns:nz,1)
        wl0=0.                  ! turn off later rewind operations
      endif

      wt=log(wl/wlbaer(ind))/log(wlbaer(3-ind)/wlbaer(ind))
      wt=max(zero,min(wt,one))
      if(debug) then
        print *,'ns,nmom,imoma ',ns,nmom,imoma
        print '(a,f8.3,a,2f8.3)','wl=',wl,'    wl1-wl2:',wlbaer(1:2)
        print '(5x,3a11)','tau','wbaer','pmom'
      endif
      do i=ns,nz
        if(min(taer(i,ind),taer(i,3-ind)).gt.0.) then
          dtau(i)=taer(i,ind)*(taer(i,3-ind)/taer(i,ind))**wt
        else
          dtau(i)=taer(i,ind)*(1.-wt)+taer(i,3-ind)*wt
        endif
        wbaer(i)=waer(i,ind)*(1.-wt)+waer(i,3-ind)*wt

        if(nmom.eq.1) then
          nmomg=maxmom
          gg=gaer(1,i,ind)*(1.-wt)+gaer(1,i,3-ind)*wt
          call getmom(imoma,gg,nmomg,pm(0))
        else
          nmomg=nmom
          pm(1:nmom)=gaer(1:nmom,i,ind)*(1.-wt)+gaer(1:nmom,i,3-ind)*wt
        endif
        pmom(1:nmomg,i)=pmom(1:nmomg,i)+pm(1:nmomg)*dtau(i)*wbaer(i)
        if(debug) print '(i5,14es11.3)',
     &       i,dtau(i),wbaer(i),pm(1:min(nmomg,12))
      enddo

      return
 20   continue
      print *,'not enough aerosol records'
      print *,'need ',nz-ns+1,' records'
      print '(a/(10es11.3))','taer',taer(ns:nz,ind)
      stop
 30   continue
      print *,'no data found in aerosol.dat'
      stop
      end
c=======================================================================
 
