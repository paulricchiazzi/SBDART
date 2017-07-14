c file:                  spectra.f
c
c external routines:     salbedo,suralb,raysig,rayleigh,zensun,solirr
c                        setfilt,filter
c
c internal routines:     rdspec,snow,clearw,lakew,seaw,sand,vegeta
c                        meteo, goese,goesw,avhr81,avhr82,avhr91,avhr92,
c                        avhr101,avhr102,avhr111,avhr112,gtr1,gtr2,nm410
c                        nm936,usersat,mfrsr1,mfrsr2,mfrsr3,mfrsr4,mfrsr5,
c                        mfrsr6,avhr83,avhr84,avhr85,setlow,sun1s,sunlow,
c                        sunmod
c
c internal modules:      albblk,  sunblk, fltblk
c=======================================================================
      module albblk             ! albedo data base
      use params, only: kr
      implicit none
      save
      integer, parameter :: mxwl=1000, mxbc=101
      real(kr), dimension(mxwl) :: wlalb,alb
      real(kr) :: wndspd,wndang,chlor,salin
      real(kr) :: hssa,hasym,hotspt,hotwdth
      real(kr) :: rliso,rlvol,rlgeo,rlhot,rlwdth

      integer :: nna=0,ibdrf
      end module albblk
c-----------------------------------------------------------------------
      function salbedo(wl)
c 
c input:   wl
c output:  surface albedo
c
      use params, only : one, zero, kr
      use albblk, only : nna, alb, wlalb
      implicit none
      integer :: j
      real(kr) :: wl, wt, salbedo
c mfl debugging
      character (len=7)  :: num1
      character (len=4)  :: num2

      if(nna.eq.0) return ! such as for isalb eq 7

c mfl debugging
c     print *, 'asdfasdfadfadfad'

      if(wl.lt.wlalb(1))
     &  call errmsg(18,'SALBEDO--spectral range lt')

c  mfl
      if(wl.lt.wlalb(1)) then
c       print *, 'WARNING: SALBEDO--spectral range lt', wlalb(1), wl
      endif
          
      write(num1, '(f7.1)') wlalb(nna)
      write(num2, '(i4)') nna

      if(wl.gt.wlalb(nna))
     &     call errmsg(18,'SALBEDO--spectral range gt '//num1)

c  mfl
      if(wl.gt.wlalb(nna)) then
c       print *, 'WARNING: SALBEDO--spectral range gt', wlalb(nna)
      endif
      
      call locate(wlalb,nna,wl,j)
      wt=(wl-wlalb(j))/(wlalb(j+1)-wlalb(j))
      wt=max(zero,min(one,wt))
      salbedo=alb(j)*(1.-wt)+alb(j+1)*wt
c mfl
c     print *, salbedo
      return
      end
c-----------------------------------------------------------------------
      subroutine suralb(isalb,albcon,sc)
c
c purpose: specify surface composition.  Must be called before SALBEDO
c
c input:   isalb    surface type
c                    -1 -spectral surface albedo read from "albedo.dat"
c                     0 -user specified, spectrally uniform, surface albedo
c                     1 -snow
c                     2 -clear water
c                     3 -lake water
c                     4 -sea  water
c                     5 -sand
c                     6 -vegetation
c                     7 -ocean bdrf model (from 6s)
c                     8 -hapke bdrf model
c                    10 -combination of snow, seawater, sand and vegetation
c
c         isalb
c           0    albcon --  constant surface albedo
c
c           7    sc     --  ocean reflection parameters 
c                           sc(1)=Oceanic pigment concentration in mg/m3
c                           sc(2)=wind speed in meters/second
c                           sc(3)=wind direction azimuth angle (degrees)
c                           sc(4)=oceanic salinity in parts per thousand
c
c           8    sc     --  hapke bdrf model parameters 
c                           sc(1)=surface particle single scattering albedo
c                           sc(2)=surface particle asymmetry factor
c                           sc(3)=hot spot magnitude
c                           sc(4)=width of hotspot
c
c           9    sc     --  Ross-thick Li-sparse bdrf model parameters 
c                           sc(1)=isotropic coefficient
c                           sc(2)=volumetric coefficient
c                           sc(3)=geometric shadowing coefficient
c                           sc(4)=hot spot magnitude
c                           sc(5)=hot spot width
c
c          10   sc      --  surface component fractions
c                           sc(1)=fraction of snow
c                           sc(2)=fraction of ocean
c                           sc(3)=fraction of sand
c                           sc(4)=fraction of vegetation
c

      use params, only: kr
      use albblk, only: mxwl, nna, alb, wlalb,ibdrf,
     &     wndspd, wndang, salin, chlor, 
     &     hssa, hasym, hotspt, hotwdth,
     &     rliso,rlvol,rlgeo,rlhot,rlwdth
      implicit none
      integer :: isalb
      real(kr) :: albcon, sc(5), albx(mxwl), mui

      nna=mxwl

c  mfl     if(isalb.eq.-1) then
c        return
c      endif

c mfl
c     print *, 'isalb = ', isalb

      select case (isalb)
      case (-1) ; call rdspec('albedo.dat',wlalb,alb,nna)
      case (0)  
        alb(1:2)=albcon
        wlalb(1:2)=(/0._kr,huge(0._kr)/)
        nna=2

      case (1)  ; call snow(wlalb,alb,nna)
      case (2)  ; call clearw(wlalb,alb,nna)
      case (3)  ; call lakew(wlalb,alb,nna)
      case (4)  ; call seaw(wlalb,alb,nna)
      case (5)  ; call sand(wlalb,alb,nna)
      case (6)  ; call vegeta(wlalb,alb,nna)

      case (-7,7)                ! ocean reflectance model
        ibdrf=1
        nna=0
        chlor=sc(1)
        wndspd=sc(2)
        salin=sc(3)

      case (-8,8)                ! hapke bdrf model
        ibdrf=2
        nna=0
        hssa=sc(1)
        hasym=sc(2)
        hotspt=sc(3)
        hotwdth=sc(4)

      case (-9,9)                ! ross-thick, li-sparse bdrf model
        ibdrf=3
        nna=0
        rliso=sc(1)
        rlvol=sc(2)
        rlgeo=sc(3)
        rlhot=sc(4)
        rlwdth=sc(5)

      case (10)               ! spectral mixture
        call snow(wlalb,albx,nna)   ; alb=albx*sc(1)
        call seaw(wlalb,albx,nna)   ; alb=albx*sc(2)+alb
        call sand(wlalb,albx,nna)   ; alb=albx*sc(3)+alb
        call vegeta(wlalb,albx,nna) ; alb=albx*sc(4)+alb

      case default 
        print *,'ERROR in suralb --- illegal value of isalb',isalb
        stop
      end select

      return
      end

c=======================================================================
      function raysig(v)

c purpose:
c    calculate molecular rayleigh scattering coefficient
c    using approximation of shettle et al., 1980 (appl. opt., 2873-4)
c    with the depolarization = 0.0279 instead of 0.035
c    for temperature = 273 k & pressure = 1 atm.
c
c input:
c  v         wavenumber cm-1
c
c output:
c  raysig    scattering coefficient (km-1) 
c            optical depth = raysig * (p/pzero)*(tzero/t)*dz

      use params, only: kr
      implicit none
      real(kr), parameter :: fit1=9.38076e+18,
     &                       fit2=-1.08426e+09
      real(kr) :: v, raysig


      raysig = v**4/(fit1+fit2*v**2)
      return
      end
c=======================================================================

      subroutine rayleigh(wl,z,p,t,nz,dtaur)

c  purpose:
c
c  input:  
c    wl        wavelength in microns
c    z         altitude array, z(1)=0 (km)
c    p         pressure array, p(1) at surface (millibars)
c    t         temperature array, t(1) at surface (kelvin)
c    nz        number of atmospheric layers
c
c  output: 
c    dtaur    increments of rayleigh scattering optical depth
c             dtaur(nz) represents the optical depth of the bottom layer
c
      use params, only: mxly, one, zero, pzero, tzero, kr
      implicit none
      integer :: i, im, nz 
      real(kr) :: z(*),p(*),t(*),dtaur(*),raysig,sig, dz, rhom, rhop, wl
c
      dtaur(1:nz)=0.

c rayleigh scattering coefficient (1/km) 

      sig=raysig(10000./wl)

      dtaur(1)=sig*(p(nz)/pzero)/(t(nz)/tzero)*5. ! assume 5km scale ht.

      do i=2,nz
        im=nz-i+1
        rhom=(p(im)/pzero)/(t(im)/tzero)
        rhop=(p(im+1)/pzero)/(t(im+1)/tzero)
        dz=z(im+1)-z(im)
        if(rhom.eq.rhop) then
          dtaur(i)=.5*sig*dz*(rhom+rhop)
        else
          dtaur(i)=sig*dz*(rhop-rhom)/log(rhop/rhom)
        endif
      enddo

      return
      end
c=======================================================================
       function bdref(wvnmlo, wvnmhi, mur, mui, phir )

c      Supplies surface bi-directional reflectivity.
c
c      NOTE 1: Bidirectional reflectivity in DISORT is defined
c              by Eq. 39 in STWL.
c      NOTE 2: Both MUR and MUI (cosines of reflection and incidence
c              angles) are positive.
c
c  INPUT:
c
c    WVNMLO : Lower wavenumber (inv cm) of spectral interval
c
c    WVNMHI : Upper wavenumber (inv cm) of spectral interval
c
c    MUR    : Cosine of angle of reflection (positive)
c
c    MUI    : Cosine of angle of incidence (positive)
c
c    PHIR   : Difference of azimuth angles of incidence and reflection
c                (radians)
c
c
c   Called by- DREF, SURFAC, BDREFCHK

c +-------------------------------------------------------------------+
c
c     .. Scalar Arguments ..

      use albblk, only: ibdrf
      use params, only: kr

      implicit none
      real(kr) :: phid, phir, mur, mui, wvnmhi, wvnmlo, wl, bdref, blim

      wl=20000./(wvnmhi+wvnmlo)
      select case (ibdrf)
      case (1); call seabdrf(wl,mui,mur,phir,bdref)
      case (2); call hapkbdrf(mui,mur,phir,bdref)
      case (3); call rtlsbdrf(mui,mur,phir,bdref)

      case default
        print *, 'Error in BDREF --- illegal value of ibdrf',ibdrf
        stop
      end select

      return
      end
c=======================================================================
      subroutine hapkbdrf(ui,ur,phir,bdrf)

c returns hapke bdrf factor 
c (http://stratus.ssec.wisc.edu/streamer_web/userman/surfalb.html) 
c
c   ui    cosine of incident angle (positive)
c   ur    cosine of reflection angle (positive)
c   phir  difference of solar and viewing azimuth (radians)
c
c output:
c   bdref surface reflectance
 
      use params, only: kr, pi
      use albblk, only:
     &     hssa,      ! ground particle single scattering albedo (0,1)
     &     hasym,     ! ground particle asymmetry factor         (-1,1)
     &     hotspt,    ! hot spot amplitude factor 
     &     hotwdth    ! hot spot width factor

      implicit none
      
      real(kr) ::
     &     bdrf,                ! returned bdrf value
     &     s,                   ! angle between incomming and reflected ray
     &     coss,                ! cosine of s
     &     ur,                  ! cosine of reflected zenith angle
     &     ui,                  ! cosine of incident zenith angle
     &     phir,                ! azimuth between incident and reflected rays
     &     pfun,                ! henyey greenstein scattering function
     &     pfun0,               ! h.g. function at coss=1
     &     b0,                  ! hot spot factor
     &     bfun,                ! hot spot function
     &     hfunr,               ! multiple scattering function (reflection)
     &     hfuni                ! multiple scattering function (incidence)

      coss=ui * ur +sqrt(1.-ur**2)*sqrt(1.-ui**2)*cos(pi-phir)

      s=acos(coss)

      pfun=(1.-hasym**2)/(1+hasym**2+2*hasym*coss)**1.5
      pfun0=(1.-hasym**2)/(1+hasym)**3
      b0=hotspt/(hssa*pfun0)
      bfun=b0/(1.+tan(s/2)/hotwdth)
      hfunr=(1.+2*ur)/(1.+2.*ur*sqrt(1.-hssa))
      hfuni=(1.+2*ui)/(1.+2.*ui*sqrt(1.-hssa))

      bdrf=(1.+bfun)*pfun + hfunr*hfuni - 1.
      bdrf=.25*hssa*bdrf/(ur+ui)

      return
      end
c=======================================================================
      subroutine rtlsbdrf(mui,mur,phir,bdrf)

c returns ross-thick, li-sparse  bdrf factor 
c (Ross, 1981; Li and Strahler, 1992).
c (http://stratus.ssec.wisc.edu/streamer_web/userman/surfalb.html) 
c
c   mui    cosine of incident angle (positive)
c   mur    cosine of reflection angle (positive)
c   phir  difference of solar and viewing azimuth (radians)
c
c output:
c   bdref surface reflectance
 
      use params, only: kr, pi
      use albblk, only:
     &     rliso,               ! isotropic coefficient
     &     rlvol,               ! volumetric coefficient
     &     rlgeo,               ! geometric coefficient
     &     rlhot,               ! hot spot magnitude
     &     rlwdth               ! hot spot width

      implicit none
      
      real(kr) :: sp, vza, sza, tanvzap, tanszap, vzap, szap, cossp,
     &     dd, secsum, cost, t, coss, sins, s, f1, f2, mui, mur,
     &     ui, ur, phir, bdrf, cosra

      ui=max(mui,.01_kr)
      ur=max(mur,.01_kr)
      cosra=cos(pi-phir)
      coss=ui * ur +sqrt(1.-ur**2)*sqrt(1.-ui**2)*cosra
      coss=max(-1._kr,min(coss,1._kr))
      s=acos(coss)
      sins=sin(s)

                                ! thick canopy reflectance kernel
      f1=(pi/2-s)*coss+sins     
      f1=f1/(ui+ur)-pi/4.

      vza=acos(ur)
      sza=acos(ui)
      tanvzap=rlwdth*tan(vza)
      tanszap=rlwdth*tan(sza)
      if(rlwdth.eq.1._kr) then
        vzap=vza
        szap=sza
      else
        vzap=atan(tanvzap)
        szap=atan(tanszap)
      endif
      
      cossp=cos(szap)*cos(vzap)+sin(szap)*sin(vzap)*cosra
      cossp=max(-1._kr,min(cossp,1._kr))
      dd=tanszap**2+tanvzap**2-2*tanszap*tanvzap*cosra

      secsum=1./cos(szap)+1./cos(vzap)
      cost=rlhot*sqrt(dd+(tanszap*tanvzap*sin(pi-phir))**2)
      cost=cost/secsum
      cost=max(-1._kr,min(cost,1._kr))

      t=acos(cost)
                                ! sparse object reflectance kernel
      
      f2=(t-sin(t)*cost)*secsum/pi ! sparse object kernel
      f2=f2-1./cos(vzap)+.5*(1.+cossp)/(cos(szap)*cos(vzap))

      bdrf=rliso+rlvol*f1+rlgeo*f2

      return
      end
c=======================================================================
      subroutine seabdrf(wl,mus,muv,phir,bdref)

c returns reflectance of sea surface due to surface reflection,
c foam and subsurface particulates and rayleigh scattering
c
c input:
c   wl    wavelength (um)
c   mus   cosine of solar zenith angle (0-1)
c   muv   cosine of viewing zenith angle (0-1)
c   phir  difference of solar and viewing azimuth (radians)
c
c output:
c   bdref surface reflectance
 
      use params, only: kr, pi
      use albblk, only:
     &     wndspd,    ! wind speed in m/s
     &     wndang,    ! solar_azimuth - wind_azimuth (degrees)
     &     chlor,     ! chlorphyll conentration (mg/m3) 
     &     salin      ! salinity (ppt) sea water 32-37  fresh water < .5
      implicit none
      
      real(kr) :: wl,ts,tv,sunaz,rgl,bdref,mus,muv,phir
      real(kr), save :: wllast=0., wndwt=-1., nr, ni, rsw, rfoam

      ! foam parameters from koepke, 1984

      real(kr), parameter ::
     &     wndc1=2.951e-6, ! coefficient 1 of foam area formula
     &     wndc2=3.52,     ! coefficient 2 of foam area formula
     &     rfco=0.22       ! average albedo of foam +- 0.11 (koepke)

      if(wllast.ne.wl) then 
        call indwat(wl,chlor,nr,ni)   ! index of refraction
        call morcasiwat(wl,chlor,rsw) ! sub-surface contribution (lambertian)
        if(chlor .eq. 0.) rsw=0.
        wndwt=wndc1*wndspd**wndc2     ! area covered by foam
        rfoam=wndwt*rfco              ! foam reflectance (lambertian)
        wllast=wl
      endif
      call sunglint(wndspd,nr,ni,mus,muv,phir,rgl)
      bdref=rfoam+(1.-wndwt)*rgl+(1.-rfoam)*rsw
      return
      end
c=======================================================================
            
      subroutine morcasiwat(wl,c,rsw)
c 
c spectral diffuse attenuation coefficient of case i waters as predicted 
c by morel within the spectral range 400-700nm (1988, journal of geophysical 
c research, vol.93, no c9, pp 10749-10768)
c
c this subroutine copied whole from tanre's 6s code.
c
c input parameters:     wl wavelength (in micrometers)
c                       c  pigment concentration 
c                          (chlorophyll + pheophytin mg/m3) 
c output parameter:     rsw  reflectance of water
c
c according morel,1988, we use:
c
c kd    spectral value of the attenuation coefficient for 
c        downwelling irradiance
c        with: kd=kw+xc*c**e
c kw    spectral value of the diffuse attenuation coefficient 
c        for pure oceanic water
c xc, e spectral coefficients to compute the diffuse attenuation 
c        coefficient for pigment
c bb    total backscattering coefficient
c        with: bb=0.5*bw+bbt*b
c bw    spectral value of the molecular scattering coefficient of water
c bbt,b parameters to compute the scattering coefficients of pigments
c
c rsw   subsurface contribution to ocean reflectivity (transfered to above
c       the surface by averaging fresnel transmission through surface)
c
c        with: rsw=(0.33/u)*(bb/kd)     where u is depends on rsw
c
      use params, only: kr
      implicit none
      real(kr) :: kw,kd,tkw(61),txc(61),te(61),tbw(61),
     &     wl,c,rsw,xc,e,bw,bb,b,bbt,u1,r1,u2,err
      integer :: iwl

      data tkw/
     & 0.0209,0.0200,0.0196,0.0189,0.0183,
     & 0.0182,0.0171,0.0170,0.0168,0.0166,
     & 0.0168,0.0170,0.0173,0.0174,0.0175,
     & 0.0184,0.0194,0.0203,0.0217,0.0240,
     & 0.0271,0.0320,0.0384,0.0445,0.0490,
     & 0.0505,0.0518,0.0543,0.0568,0.0615,
     & 0.0640,0.0640,0.0717,0.0762,0.0807,
     & 0.0940,0.1070,0.1280,0.1570,0.2000,
     & 0.2530,0.2790,0.2960,0.3030,0.3100,
     & 0.3150,0.3200,0.3250,0.3300,0.3400,
     & 0.3500,0.3700,0.4050,0.4180,0.4300,
     & 0.4400,0.4500,0.4700,0.5000,0.5500,
     & 0.6500/
      data txc/
     & 0.1100,0.1110,0.1125,0.1135,0.1126,
     & 0.1104,0.1078,0.1065,0.1041,0.0996,
     & 0.0971,0.0939,0.0896,0.0859,0.0823,
     & 0.0788,0.0746,0.0726,0.0690,0.0660,
     & 0.0636,0.0600,0.0578,0.0540,0.0498,
     & 0.0475,0.0467,0.0450,0.0440,0.0426,
     & 0.0410,0.0400,0.0390,0.0375,0.0360,
     & 0.0340,0.0330,0.0328,0.0325,0.0330,
     & 0.0340,0.0350,0.0360,0.0375,0.0385,
     & 0.0400,0.0420,0.0430,0.0440,0.0445,
     & 0.0450,0.0460,0.0475,0.0490,0.0515,
     & 0.0520,0.0505,0.0440,0.0390,0.0340,
     & 0.0300/
      data te/
     & 0.668,0.672,0.680,0.687,0.693,
     & 0.701,0.707,0.708,0.707,0.704,
     & 0.701,0.699,0.700,0.703,0.703,
     & 0.703,0.703,0.704,0.702,0.700,
     & 0.700,0.695,0.690,0.685,0.680,
     & 0.675,0.670,0.665,0.660,0.655,
     & 0.650,0.645,0.640,0.630,0.623,
     & 0.615,0.610,0.614,0.618,0.622,
     & 0.626,0.630,0.634,0.638,0.642,
     & 0.647,0.653,0.658,0.663,0.667,
     & 0.672,0.677,0.682,0.687,0.695,
     & 0.697,0.693,0.665,0.640,0.620,
     & 0.600/
      data tbw/
     & 0.0076,0.0072,0.0068,0.0064,0.0061,
     & 0.0058,0.0055,0.0052,0.0049,0.0047,
     & 0.0045,0.0043,0.0041,0.0039,0.0037,
     & 0.0036,0.0034,0.0033,0.0031,0.0030,
     & 0.0029,0.0027,0.0026,0.0025,0.0024,
     & 0.0023,0.0022,0.0022,0.0021,0.0020,
     & 0.0019,0.0018,0.0018,0.0017,0.0017,
     & 0.0016,0.0016,0.0015,0.0015,0.0014,
     & 0.0014,0.0013,0.0013,0.0012,0.0012,
     & 0.0011,0.0011,0.0010,0.0010,0.0010,
     & 0.0010,0.0009,0.0008,0.0008,0.0008,
     & 0.0007,0.0007,0.0007,0.0007,0.0007,
     & 0.0007/
      if (wl.lt.0.400.or.wl.gt.0.700)then
        rsw=0.000
        goto 60
      endif

      iwl=1+nint((wl-0.400)/0.005)
      kw=tkw(iwl)
      xc=txc(iwl)
      e=te(iwl)
      bw=tbw(iwl)
c
      if (abs(c).lt.0.0001)then
        bb=0.5*bw
        kd=kw
      else
        b=0.30*c**0.62
        bbt=0.002+0.02*(0.5-0.25*log10(c))*0.550/wl
        bb=0.5*bw+bbt*b
        kd=kw+xc*c**e
      endif

      u1=0.75
      r1=0.33*bb/u1/kd

 50   u2=0.90*(1.-r1)/(1.+2.25*r1)
      rsw=0.33*bb/u2/kd
      err=abs((rsw-r1)/rsw)
      if (err.lt.0.0001)goto 60
      r1=rsw
      goto 50
 60   return
      end
c=======================================================================
       subroutine indwat(wl,xsal,nr,ni)
c
c input parameters:  wl=wavelength (in micrometers)
c                    xsal=salinity (in ppt), if xsal<0 then 34.3ppt by default
c output parameters: nr=index of refraction of sea water
c                    ni=extinction coefficient of sea water
c
c temperature dependence of ni ignored
c 
c
       use params, only: kr
       implicit none

       integer, parameter :: ntab=860
       integer :: i
       real(kr), parameter :: nrc=0.006, nic=0.000
       real(kr) ::  wltab(ntab),mrtab(ntab),mitab(ntab),nr,ni,wl,xsal,wt

c      (1) for 0.01 - 1.e7 microns :  Segelstein, D., 1981:
c             "The Complex Refractive Index of Water", M.S. Thesis,
c             University of Missouri--Kansas City

c      (2) for 10. - 1.e7 microns:  Ray, P., 1972:  Broadband Complex
c             Refractive Indices of Ice and Water, Appl. Opt. 11,
c             1836-1844

c      (There is a new reference, Wieliczka, D. et al., Appl. Opt.
c       28, 1714-1719, 1989, with some updated data for the IR)

      data wltab(1:50)/
     &  1.000e-02, 1.099e-02, 1.200e-02, 1.300e-02, 1.400e-02,
     &  1.600e-02, 1.799e-02, 2.000e-02, 2.198e-02, 2.399e-02,
     &  2.600e-02, 2.799e-02, 2.999e-02, 3.199e-02, 3.396e-02,
     &  3.598e-02, 3.802e-02, 3.999e-02, 4.198e-02, 4.395e-02,
     &  4.603e-02, 4.797e-02, 5.000e-02, 5.200e-02, 5.395e-02,
     &  5.598e-02, 5.794e-02, 5.998e-02, 6.194e-02, 6.397e-02,
     &  6.607e-02, 6.808e-02, 6.998e-02, 7.194e-02, 7.396e-02,
     &  7.603e-02, 7.798e-02, 7.998e-02, 8.203e-02, 8.395e-02,
     &  8.590e-02, 8.790e-02, 8.995e-02, 9.205e-02, 9.397e-02,
     &  9.594e-02, 9.795e-02, 1.000e-01, 1.021e-01, 1.040e-01/

      data wltab(51:100)/
     &  1.059e-01, 1.079e-01, 1.099e-01, 1.119e-01, 1.140e-01,
     &  1.159e-01, 1.180e-01, 1.200e-01, 1.219e-01, 1.239e-01,
     &  1.259e-01, 1.279e-01, 1.300e-01, 1.321e-01, 1.340e-01,
     &  1.361e-01, 1.380e-01, 1.400e-01, 1.419e-01, 1.439e-01,
     &  1.459e-01, 1.479e-01, 1.500e-01, 1.520e-01, 1.542e-01,
     &  1.560e-01, 1.581e-01, 1.600e-01, 1.622e-01, 1.641e-01,
     &  1.660e-01, 1.679e-01, 1.698e-01, 1.722e-01, 1.742e-01,
     &  1.750e-01, 1.799e-01, 1.849e-01, 1.901e-01, 1.950e-01,
     &  2.000e-01, 2.051e-01, 2.099e-01, 2.148e-01, 2.198e-01,
     &  2.249e-01, 2.301e-01, 2.350e-01, 2.399e-01, 2.449e-01/
      data wltab(101:150)/
     &  2.500e-01, 2.553e-01, 2.600e-01, 2.648e-01, 2.698e-01,
     &  2.748e-01, 2.799e-01, 2.851e-01, 2.897e-01, 2.951e-01,
     &  2.999e-01, 3.048e-01, 3.097e-01, 3.148e-01, 3.199e-01,
     &  3.251e-01, 3.304e-01, 3.350e-01, 3.396e-01, 3.451e-01,
     &  3.499e-01, 3.548e-01, 3.598e-01, 3.648e-01, 3.698e-01,
     &  3.750e-01, 3.802e-01, 3.846e-01, 3.899e-01, 3.954e-01,
     &  3.999e-01, 4.046e-01, 4.102e-01, 4.149e-01, 4.198e-01,
     &  4.246e-01, 4.295e-01, 4.345e-01, 4.395e-01, 4.446e-01,
     &  4.498e-01, 4.550e-01, 4.603e-01, 4.645e-01, 4.699e-01,
     &  4.753e-01, 4.797e-01, 4.853e-01, 4.898e-01, 4.955e-01/
      data wltab(151:200)/
     &  5.000e-01, 5.047e-01, 5.105e-01, 5.152e-01, 5.200e-01,
     &  5.248e-01, 5.297e-01, 5.346e-01, 5.395e-01, 5.445e-01,
     &  5.495e-01, 5.546e-01, 5.598e-01, 5.649e-01, 5.702e-01,
     &  5.754e-01, 5.794e-01, 5.848e-01, 5.902e-01, 5.957e-01,
     &  5.998e-01, 6.053e-01, 6.095e-01, 6.152e-01, 6.194e-01,
     &  6.252e-01, 6.295e-01, 6.353e-01, 6.397e-01, 6.456e-01,
     &  6.501e-01, 6.546e-01, 6.607e-01, 6.653e-01, 6.699e-01,
     &  6.745e-01, 6.808e-01, 6.855e-01, 6.902e-01, 6.950e-01,
     &  6.998e-01, 7.047e-01, 7.096e-01, 7.145e-01, 7.195e-01,
     &  7.244e-01, 7.295e-01, 7.345e-01, 7.396e-01, 7.447e-01/
      data wltab(201:250)/
     &  7.499e-01, 7.551e-01, 7.603e-01, 7.656e-01, 7.691e-01,
     &  7.745e-01, 7.798e-01, 7.852e-01, 7.907e-01, 7.943e-01,
     &  7.998e-01, 8.054e-01, 8.091e-01, 8.147e-01, 8.204e-01,
     &  8.241e-01, 8.298e-01, 8.356e-01, 8.395e-01, 8.453e-01,
     &  8.492e-01, 8.551e-01, 8.590e-01, 8.650e-01, 8.710e-01,
     &  8.750e-01, 8.790e-01, 8.851e-01, 8.892e-01, 8.954e-01,
     &  8.995e-01, 9.057e-01, 9.099e-01, 9.141e-01, 9.204e-01,
     &  9.247e-01, 9.290e-01, 9.354e-01, 9.397e-01, 9.441e-01,
     &  9.506e-01, 9.550e-01, 9.594e-01, 9.660e-01, 9.705e-01,
     &  9.750e-01, 9.795e-01, 9.840e-01, 9.908e-01, 9.954e-01/
      data wltab(251:300)/
     &  1.000e+00, 1.009e+00, 1.021e+00, 1.030e+00, 1.040e+00,
     &  1.050e+00, 1.059e+00, 1.069e+00, 1.079e+00, 1.089e+00,
     &  1.099e+00, 1.109e+00, 1.119e+00, 1.130e+00, 1.140e+00,
     &  1.151e+00, 1.159e+00, 1.169e+00, 1.180e+00, 1.191e+00,
     &  1.200e+00, 1.211e+00, 1.219e+00, 1.230e+00, 1.239e+00,
     &  1.250e+00, 1.259e+00, 1.271e+00, 1.279e+00, 1.291e+00,
     &  1.300e+00, 1.309e+00, 1.321e+00, 1.331e+00, 1.340e+00,
     &  1.349e+00, 1.361e+00, 1.371e+00, 1.380e+00, 1.390e+00,
     &  1.400e+00, 1.409e+00, 1.419e+00, 1.429e+00, 1.439e+00,
     &  1.449e+00, 1.459e+00, 1.469e+00, 1.479e+00, 1.489e+00/
      data wltab(301:350)/
     &  1.500e+00, 1.510e+00, 1.520e+00, 1.531e+00, 1.542e+00,
     &  1.549e+00, 1.560e+00, 1.570e+00, 1.581e+00, 1.589e+00,
     &  1.600e+00, 1.611e+00, 1.622e+00, 1.629e+00, 1.641e+00,
     &  1.648e+00, 1.660e+00, 1.671e+00, 1.679e+00, 1.690e+00,
     &  1.698e+00, 1.710e+00, 1.722e+00, 1.730e+00, 1.742e+00,
     &  1.750e+00, 1.762e+00, 1.770e+00, 1.778e+00, 1.791e+00,
     &  1.799e+00, 1.811e+00, 1.820e+00, 1.828e+00, 1.841e+00,
     &  1.849e+00, 1.862e+00, 1.871e+00, 1.879e+00, 1.888e+00,
     &  1.901e+00, 1.910e+00, 1.919e+00, 1.932e+00, 1.941e+00,
     &  1.950e+00, 1.959e+00, 1.968e+00, 1.982e+00, 1.991e+00/
      data wltab(351:400)/
     &  2.000e+00, 2.009e+00, 2.018e+00, 2.028e+00, 2.042e+00,
     &  2.051e+00, 2.061e+00, 2.070e+00, 2.080e+00, 2.089e+00,
     &  2.099e+00, 2.109e+00, 2.118e+00, 2.128e+00, 2.138e+00,
     &  2.148e+00, 2.158e+00, 2.168e+00, 2.178e+00, 2.188e+00,
     &  2.198e+00, 2.208e+00, 2.218e+00, 2.228e+00, 2.239e+00,
     &  2.249e+00, 2.259e+00, 2.270e+00, 2.280e+00, 2.291e+00,
     &  2.301e+00, 2.312e+00, 2.317e+00, 2.328e+00, 2.339e+00,
     &  2.350e+00, 2.361e+00, 2.371e+00, 2.382e+00, 2.388e+00,
     &  2.399e+00, 2.410e+00, 2.421e+00, 2.432e+00, 2.438e+00,
     &  2.449e+00, 2.460e+00, 2.472e+00, 2.477e+00, 2.489e+00/
      data wltab(401:450)/
     &  2.500e+00, 2.512e+00, 2.518e+00, 2.529e+00, 2.541e+00,
     &  2.553e+00, 2.564e+00, 2.570e+00, 2.576e+00, 2.582e+00,
     &  2.588e+00, 2.594e+00, 2.606e+00, 2.612e+00, 2.618e+00,
     &  2.624e+00, 2.630e+00, 2.636e+00, 2.648e+00, 2.655e+00,
     &  2.661e+00, 2.667e+00, 2.673e+00, 2.679e+00, 2.685e+00,
     &  2.698e+00, 2.704e+00, 2.710e+00, 2.716e+00, 2.723e+00,
     &  2.729e+00, 2.742e+00, 2.748e+00, 2.754e+00, 2.761e+00,
     &  2.767e+00, 2.780e+00, 2.786e+00, 2.792e+00, 2.799e+00,
     &  2.812e+00, 2.818e+00, 2.825e+00, 2.831e+00, 2.838e+00,
     &  2.851e+00, 2.858e+00, 2.864e+00, 2.871e+00, 2.884e+00/
      data wltab(451:500)/
     &  2.891e+00, 2.897e+00, 2.904e+00, 2.917e+00, 2.924e+00,
     &  2.931e+00, 2.938e+00, 2.951e+00, 2.958e+00, 2.965e+00,
     &  2.978e+00, 2.985e+00, 2.999e+00, 3.048e+00, 3.097e+00,
     &  3.148e+00, 3.199e+00, 3.251e+00, 3.304e+00, 3.350e+00,
     &  3.396e+00, 3.451e+00, 3.499e+00, 3.548e+00, 3.598e+00,
     &  3.647e+00, 3.698e+00, 3.750e+00, 3.802e+00, 3.846e+00,
     &  3.899e+00, 3.954e+00, 3.999e+00, 4.046e+00, 4.102e+00,
     &  4.149e+00, 4.198e+00, 4.246e+00, 4.295e+00, 4.345e+00,
     &  4.395e+00, 4.446e+00, 4.498e+00, 4.550e+00, 4.603e+00,
     &  4.645e+00, 4.699e+00, 4.753e+00, 4.797e+00, 4.853e+00/
      data wltab(501:550)/
     &  4.898e+00, 4.955e+00, 5.000e+00, 5.047e+00, 5.105e+00,
     &  5.152e+00, 5.200e+00, 5.248e+00, 5.297e+00, 5.346e+00,
     &  5.395e+00, 5.445e+00, 5.495e+00, 5.546e+00, 5.598e+00,
     &  5.649e+00, 5.702e+00, 5.754e+00, 5.794e+00, 5.848e+00,
     &  5.902e+00, 5.957e+00, 5.998e+00, 6.053e+00, 6.095e+00,
     &  6.152e+00, 6.194e+00, 6.252e+00, 6.295e+00, 6.353e+00,
     &  6.397e+00, 6.457e+00, 6.501e+00, 6.546e+00, 6.607e+00,
     &  6.653e+00, 6.699e+00, 6.745e+00, 6.808e+00, 6.855e+00,
     &  6.902e+00, 6.950e+00, 6.998e+00, 7.047e+00, 7.096e+00,
     &  7.145e+00, 7.194e+00, 7.244e+00, 7.295e+00, 7.345e+00/
      data wltab(551:600)/
     &  7.396e+00, 7.447e+00, 7.499e+00, 7.551e+00, 7.603e+00,
     &  7.656e+00, 7.691e+00, 7.745e+00, 7.798e+00, 7.852e+00,
     &  7.907e+00, 7.943e+00, 7.998e+00, 8.054e+00, 8.091e+00,
     &  8.147e+00, 8.204e+00, 8.241e+00, 8.299e+00, 8.356e+00,
     &  8.395e+00, 8.453e+00, 8.492e+00, 8.551e+00, 8.590e+00,
     &  8.650e+00, 8.710e+00, 8.750e+00, 8.790e+00, 8.851e+00,
     &  8.892e+00, 8.954e+00, 8.995e+00, 9.057e+00, 9.099e+00,
     &  9.141e+00, 9.204e+00, 9.247e+00, 9.290e+00, 9.354e+00,
     &  9.397e+00, 9.441e+00, 9.506e+00, 9.550e+00, 9.594e+00,
     &  9.661e+00, 9.705e+00, 9.750e+00, 9.795e+00, 9.840e+00/
      data wltab(601:650)/
     &  9.908e+00, 9.954e+00, 1.000e+01, 1.005e+01, 1.009e+01,
     &  1.014e+01, 1.021e+01, 1.026e+01, 1.030e+01, 1.035e+01,
     &  1.040e+01, 1.045e+01, 1.049e+01, 1.054e+01, 1.059e+01,
     &  1.064e+01, 1.069e+01, 1.074e+01, 1.079e+01, 1.084e+01,
     &  1.089e+01, 1.094e+01, 1.099e+01, 1.104e+01, 1.109e+01,
     &  1.114e+01, 1.119e+01, 1.125e+01, 1.130e+01, 1.135e+01,
     &  1.140e+01, 1.146e+01, 1.151e+01, 1.156e+01, 1.159e+01,
     &  1.164e+01, 1.170e+01, 1.175e+01, 1.180e+01, 1.186e+01,
     &  1.191e+01, 1.194e+01, 1.199e+01, 1.205e+01, 1.211e+01,
     &  1.216e+01, 1.219e+01, 1.225e+01, 1.230e+01, 1.236e+01/
      data wltab(651:700)/
     &  1.239e+01, 1.245e+01, 1.250e+01, 1.256e+01, 1.259e+01,
     &  1.265e+01, 1.271e+01, 1.276e+01, 1.279e+01, 1.285e+01,
     &  1.291e+01, 1.294e+01, 1.300e+01, 1.306e+01, 1.309e+01,
     &  1.315e+01, 1.321e+01, 1.324e+01, 1.330e+01, 1.334e+01,
     &  1.340e+01, 1.346e+01, 1.349e+01, 1.355e+01, 1.361e+01,
     &  1.365e+01, 1.371e+01, 1.374e+01, 1.380e+01, 1.384e+01,
     &  1.390e+01, 1.396e+01, 1.400e+01, 1.406e+01, 1.409e+01,
     &  1.416e+01, 1.419e+01, 1.426e+01, 1.429e+01, 1.436e+01,
     &  1.439e+01, 1.445e+01, 1.449e+01, 1.455e+01, 1.459e+01,
     &  1.465e+01, 1.469e+01, 1.476e+01, 1.479e+01, 1.486e+01/
      data wltab(701:750)/
     &  1.489e+01, 1.496e+01, 1.500e+01, 1.507e+01, 1.510e+01,
     &  1.514e+01, 1.521e+01, 1.524e+01, 1.531e+01, 1.535e+01,
     &  1.542e+01, 1.545e+01, 1.549e+01, 1.556e+01, 1.560e+01,
     &  1.567e+01, 1.570e+01, 1.574e+01, 1.581e+01, 1.585e+01,
     &  1.588e+01, 1.596e+01, 1.600e+01, 1.603e+01, 1.611e+01,
     &  1.614e+01, 1.622e+01, 1.625e+01, 1.629e+01, 1.637e+01,
     &  1.641e+01, 1.644e+01, 1.648e+01, 1.656e+01, 1.660e+01,
     &  1.663e+01, 1.671e+01, 1.675e+01, 1.679e+01, 1.687e+01,
     &  1.690e+01, 1.694e+01, 1.698e+01, 1.706e+01, 1.710e+01,
     &  1.714e+01, 1.722e+01, 1.726e+01, 1.730e+01, 1.734e+01/
      data wltab(751:800)/
     &  1.742e+01, 1.746e+01, 1.750e+01, 1.754e+01, 1.762e+01,
     &  1.766e+01, 1.770e+01, 1.774e+01, 1.778e+01, 1.786e+01,
     &  1.791e+01, 1.795e+01, 1.799e+01, 1.803e+01, 1.811e+01,
     &  1.816e+01, 1.820e+01, 1.824e+01, 1.828e+01, 1.836e+01,
     &  1.841e+01, 1.845e+01, 1.849e+01, 1.854e+01, 1.862e+01,
     &  1.866e+01, 1.871e+01, 1.875e+01, 1.879e+01, 1.884e+01,
     &  1.888e+01, 1.897e+01, 1.901e+01, 1.905e+01, 1.910e+01,
     &  1.914e+01, 1.919e+01, 1.923e+01, 1.932e+01, 1.936e+01,
     &  1.941e+01, 1.945e+01, 1.950e+01, 1.954e+01, 1.959e+01,
     &  1.963e+01, 1.968e+01, 1.977e+01, 1.982e+01, 1.986e+01/
      data wltab(801:860)/
     &  1.991e+01, 1.995e+01, 2.000e+01, 2.099e+01, 2.198e+01,
     &  2.301e+01, 2.399e+01, 2.500e+01, 2.600e+01, 2.698e+01,
     &  2.799e+01, 2.897e+01, 2.999e+01, 3.097e+01, 3.199e+01,
     &  3.304e+01, 3.396e+01, 3.499e+01, 3.598e+01, 3.698e+01,
     &  3.802e+01, 3.899e+01, 3.999e+01, 4.102e+01, 4.198e+01,
     &  4.295e+01, 4.395e+01, 4.498e+01, 4.603e+01, 4.699e+01,
     &  4.797e+01, 4.898e+01, 5.000e+01, 5.200e+01, 5.395e+01,
     &  5.598e+01, 5.794e+01, 5.998e+01, 6.194e+01, 6.397e+01,
     &  6.607e+01, 6.808e+01, 6.998e+01, 7.194e+01, 7.396e+01,
     &  7.603e+01, 7.798e+01, 7.998e+01, 8.203e+01, 8.395e+01,
     &  8.590e+01, 8.790e+01, 8.995e+01, 9.205e+01, 9.397e+01,
     &  9.594e+01, 9.795e+01, 1.000e+02, 1.099e+02, 1.200e+02/

      data mrtab(1:50)/
     &  9.684e-01, 9.648e-01, 9.610e-01, 9.570e-01, 9.528e-01,
     &  9.441e-01, 9.347e-01, 9.246e-01, 9.140e-01, 9.027e-01,
     &  8.908e-01, 8.788e-01, 8.665e-01, 8.541e-01, 8.422e-01,
     &  8.307e-01, 8.198e-01, 8.100e-01, 8.023e-01, 7.977e-01,
     &  7.970e-01, 8.056e-01, 8.207e-01, 8.310e-01, 8.352e-01,
     &  8.353e-01, 8.316e-01, 8.309e-01, 8.406e-01, 8.670e-01,
     &  9.035e-01, 9.418e-01, 9.817e-01, 1.021e+00, 1.050e+00,
     &  1.069e+00, 1.088e+00, 1.112e+00, 1.141e+00, 1.173e+00,
     &  1.215e+00, 1.259e+00, 1.303e+00, 1.347e+00, 1.388e+00,
     &  1.425e+00, 1.456e+00, 1.477e+00, 1.493e+00, 1.507e+00/
      data mrtab(51:100)/
     &  1.516e+00, 1.524e+00, 1.529e+00, 1.535e+00, 1.543e+00,
     &  1.548e+00, 1.553e+00, 1.561e+00, 1.570e+00, 1.585e+00,
     &  1.606e+00, 1.627e+00, 1.634e+00, 1.619e+00, 1.586e+00,
     &  1.536e+00, 1.496e+00, 1.471e+00, 1.461e+00, 1.461e+00,
     &  1.469e+00, 1.490e+00, 1.521e+00, 1.560e+00, 1.597e+00,
     &  1.620e+00, 1.641e+00, 1.650e+00, 1.653e+00, 1.653e+00,
     &  1.647e+00, 1.635e+00, 1.606e+00, 1.568e+00, 1.549e+00,
     &  1.543e+00, 1.513e+00, 1.492e+00, 1.475e+00, 1.463e+00,
     &  1.452e+00, 1.442e+00, 1.435e+00, 1.428e+00, 1.422e+00,
     &  1.416e+00, 1.411e+00, 1.406e+00, 1.402e+00, 1.399e+00/
      data mrtab(101:150)/
     &  1.395e+00, 1.392e+00, 1.389e+00, 1.386e+00, 1.384e+00,
     &  1.381e+00, 1.379e+00, 1.377e+00, 1.375e+00, 1.373e+00,
     &  1.371e+00, 1.370e+00, 1.368e+00, 1.367e+00, 1.365e+00,
     &  1.364e+00, 1.363e+00, 1.362e+00, 1.360e+00, 1.359e+00,
     &  1.358e+00, 1.357e+00, 1.356e+00, 1.355e+00, 1.354e+00,
     &  1.354e+00, 1.353e+00, 1.352e+00, 1.351e+00, 1.350e+00,
     &  1.350e+00, 1.349e+00, 1.348e+00, 1.348e+00, 1.347e+00,
     &  1.347e+00, 1.346e+00, 1.346e+00, 1.345e+00, 1.344e+00,
     &  1.344e+00, 1.343e+00, 1.343e+00, 1.342e+00, 1.342e+00,
     &  1.341e+00, 1.341e+00, 1.341e+00, 1.340e+00, 1.340e+00/
      data mrtab(151:200)/
     &  1.339e+00, 1.339e+00, 1.339e+00, 1.338e+00, 1.338e+00,
     &  1.338e+00, 1.337e+00, 1.337e+00, 1.337e+00, 1.336e+00,
     &  1.336e+00, 1.336e+00, 1.335e+00, 1.335e+00, 1.335e+00,
     &  1.334e+00, 1.334e+00, 1.334e+00, 1.334e+00, 1.333e+00,
     &  1.333e+00, 1.333e+00, 1.333e+00, 1.332e+00, 1.332e+00,
     &  1.332e+00, 1.332e+00, 1.331e+00, 1.331e+00, 1.331e+00,
     &  1.331e+00, 1.330e+00, 1.330e+00, 1.330e+00, 1.330e+00,
     &  1.330e+00, 1.329e+00, 1.329e+00, 1.329e+00, 1.329e+00,
     &  1.329e+00, 1.329e+00, 1.328e+00, 1.328e+00, 1.328e+00,
     &  1.328e+00, 1.328e+00, 1.328e+00, 1.328e+00, 1.327e+00/
      data mrtab(201:250)/
     &  1.327e+00, 1.327e+00, 1.327e+00, 1.327e+00, 1.327e+00,
     &  1.327e+00, 1.326e+00, 1.326e+00, 1.326e+00, 1.326e+00,
     &  1.326e+00, 1.326e+00, 1.326e+00, 1.326e+00, 1.325e+00,
     &  1.325e+00, 1.325e+00, 1.325e+00, 1.325e+00, 1.325e+00,
     &  1.325e+00, 1.325e+00, 1.325e+00, 1.324e+00, 1.324e+00,
     &  1.324e+00, 1.324e+00, 1.324e+00, 1.324e+00, 1.324e+00,
     &  1.324e+00, 1.324e+00, 1.323e+00, 1.323e+00, 1.323e+00,
     &  1.323e+00, 1.323e+00, 1.323e+00, 1.323e+00, 1.323e+00,
     &  1.323e+00, 1.323e+00, 1.322e+00, 1.322e+00, 1.322e+00,
     &  1.322e+00, 1.322e+00, 1.322e+00, 1.322e+00, 1.322e+00/
      data mrtab(251:300)/
     &  1.322e+00, 1.322e+00, 1.321e+00, 1.321e+00, 1.321e+00,
     &  1.321e+00, 1.321e+00, 1.320e+00, 1.320e+00, 1.320e+00,
     &  1.320e+00, 1.320e+00, 1.319e+00, 1.319e+00, 1.319e+00,
     &  1.319e+00, 1.319e+00, 1.319e+00, 1.318e+00, 1.318e+00,
     &  1.318e+00, 1.318e+00, 1.318e+00, 1.317e+00, 1.317e+00,
     &  1.317e+00, 1.317e+00, 1.317e+00, 1.316e+00, 1.316e+00,
     &  1.316e+00, 1.316e+00, 1.316e+00, 1.315e+00, 1.315e+00,
     &  1.315e+00, 1.315e+00, 1.315e+00, 1.314e+00, 1.314e+00,
     &  1.314e+00, 1.314e+00, 1.314e+00, 1.313e+00, 1.313e+00,
     &  1.313e+00, 1.313e+00, 1.313e+00, 1.313e+00, 1.312e+00/
      data mrtab(301:350)/
     &  1.312e+00, 1.312e+00, 1.312e+00, 1.311e+00, 1.311e+00,
     &  1.311e+00, 1.311e+00, 1.310e+00, 1.310e+00, 1.310e+00,
     &  1.310e+00, 1.309e+00, 1.309e+00, 1.309e+00, 1.309e+00,
     &  1.308e+00, 1.308e+00, 1.308e+00, 1.307e+00, 1.307e+00,
     &  1.307e+00, 1.306e+00, 1.306e+00, 1.306e+00, 1.305e+00,
     &  1.305e+00, 1.305e+00, 1.304e+00, 1.304e+00, 1.304e+00,
     &  1.303e+00, 1.303e+00, 1.303e+00, 1.302e+00, 1.302e+00,
     &  1.301e+00, 1.301e+00, 1.300e+00, 1.300e+00, 1.300e+00,
     &  1.299e+00, 1.299e+00, 1.299e+00, 1.299e+00, 1.299e+00,
     &  1.298e+00, 1.298e+00, 1.298e+00, 1.298e+00, 1.297e+00/
      data mrtab(351:400)/
     &  1.297e+00, 1.296e+00, 1.296e+00, 1.296e+00, 1.295e+00,
     &  1.294e+00, 1.294e+00, 1.293e+00, 1.293e+00, 1.292e+00,
     &  1.292e+00, 1.291e+00, 1.291e+00, 1.290e+00, 1.290e+00,
     &  1.289e+00, 1.288e+00, 1.288e+00, 1.287e+00, 1.286e+00,
     &  1.286e+00, 1.285e+00, 1.284e+00, 1.284e+00, 1.283e+00,
     &  1.282e+00, 1.281e+00, 1.280e+00, 1.280e+00, 1.279e+00,
     &  1.278e+00, 1.277e+00, 1.276e+00, 1.275e+00, 1.274e+00,
     &  1.273e+00, 1.272e+00, 1.271e+00, 1.270e+00, 1.269e+00,
     &  1.268e+00, 1.266e+00, 1.265e+00, 1.264e+00, 1.263e+00,
     &  1.261e+00, 1.260e+00, 1.258e+00, 1.257e+00, 1.255e+00/
      data mrtab(401:450)/
     &  1.253e+00, 1.251e+00, 1.250e+00, 1.248e+00, 1.246e+00,
     &  1.243e+00, 1.240e+00, 1.239e+00, 1.237e+00, 1.235e+00,
     &  1.234e+00, 1.232e+00, 1.228e+00, 1.225e+00, 1.223e+00,
     &  1.221e+00, 1.218e+00, 1.216e+00, 1.210e+00, 1.207e+00,
     &  1.203e+00, 1.199e+00, 1.195e+00, 1.191e+00, 1.188e+00,
     &  1.180e+00, 1.175e+00, 1.169e+00, 1.161e+00, 1.153e+00,
     &  1.150e+00, 1.142e+00, 1.136e+00, 1.133e+00, 1.132e+00,
     &  1.133e+00, 1.131e+00, 1.128e+00, 1.128e+00, 1.129e+00,
     &  1.128e+00, 1.126e+00, 1.125e+00, 1.128e+00, 1.133e+00,
     &  1.142e+00, 1.146e+00, 1.152e+00, 1.162e+00, 1.178e+00/
      data mrtab(451:500)/
     &  1.185e+00, 1.196e+00, 1.208e+00, 1.230e+00, 1.240e+00,
     &  1.252e+00, 1.264e+00, 1.286e+00, 1.298e+00, 1.308e+00,
     &  1.326e+00, 1.335e+00, 1.353e+00, 1.412e+00, 1.452e+00,
     &  1.467e+00, 1.462e+00, 1.449e+00, 1.433e+00, 1.417e+00,
     &  1.405e+00, 1.393e+00, 1.384e+00, 1.376e+00, 1.369e+00,
     &  1.363e+00, 1.357e+00, 1.352e+00, 1.347e+00, 1.344e+00,
     &  1.340e+00, 1.337e+00, 1.334e+00, 1.331e+00, 1.329e+00,
     &  1.326e+00, 1.324e+00, 1.322e+00, 1.320e+00, 1.318e+00,
     &  1.316e+00, 1.315e+00, 1.314e+00, 1.312e+00, 1.312e+00,
     &  1.312e+00, 1.311e+00, 1.311e+00, 1.311e+00, 1.310e+00/
      data mrtab(501:550)/
     &  1.309e+00, 1.307e+00, 1.306e+00, 1.304e+00, 1.302e+00,
     &  1.300e+00, 1.298e+00, 1.295e+00, 1.292e+00, 1.289e+00,
     &  1.286e+00, 1.282e+00, 1.278e+00, 1.274e+00, 1.269e+00,
     &  1.263e+00, 1.257e+00, 1.248e+00, 1.242e+00, 1.235e+00,
     &  1.229e+00, 1.232e+00, 1.243e+00, 1.268e+00, 1.295e+00,
     &  1.330e+00, 1.342e+00, 1.340e+00, 1.336e+00, 1.329e+00,
     &  1.325e+00, 1.320e+00, 1.318e+00, 1.315e+00, 1.311e+00,
     &  1.309e+00, 1.307e+00, 1.305e+00, 1.302e+00, 1.300e+00,
     &  1.298e+00, 1.297e+00, 1.295e+00, 1.294e+00, 1.292e+00,
     &  1.291e+00, 1.289e+00, 1.288e+00, 1.287e+00, 1.285e+00/
      data mrtab(551:600)/
     &  1.284e+00, 1.283e+00, 1.281e+00, 1.280e+00, 1.279e+00,
     &  1.277e+00, 1.276e+00, 1.275e+00, 1.273e+00, 1.272e+00,
     &  1.271e+00, 1.270e+00, 1.268e+00, 1.267e+00, 1.266e+00,
     &  1.264e+00, 1.263e+00, 1.261e+00, 1.260e+00, 1.258e+00,
     &  1.257e+00, 1.255e+00, 1.254e+00, 1.252e+00, 1.251e+00,
     &  1.249e+00, 1.247e+00, 1.246e+00, 1.245e+00, 1.243e+00,
     &  1.241e+00, 1.239e+00, 1.238e+00, 1.236e+00, 1.234e+00,
     &  1.233e+00, 1.230e+00, 1.229e+00, 1.227e+00, 1.224e+00,
     &  1.223e+00, 1.221e+00, 1.218e+00, 1.216e+00, 1.214e+00,
     &  1.211e+00, 1.209e+00, 1.207e+00, 1.204e+00, 1.202e+00/
      data mrtab(601:650)/
     &  1.199e+00, 1.196e+00, 1.193e+00, 1.190e+00, 1.187e+00,
     &  1.184e+00, 1.181e+00, 1.178e+00, 1.174e+00, 1.171e+00,
     &  1.167e+00, 1.164e+00, 1.161e+00, 1.157e+00, 1.154e+00,
     &  1.150e+00, 1.147e+00, 1.144e+00, 1.140e+00, 1.137e+00,
     &  1.134e+00, 1.131e+00, 1.129e+00, 1.125e+00, 1.122e+00,
     &  1.119e+00, 1.116e+00, 1.113e+00, 1.110e+00, 1.108e+00,
     &  1.105e+00, 1.103e+00, 1.101e+00, 1.098e+00, 1.097e+00,
     &  1.096e+00, 1.094e+00, 1.092e+00, 1.091e+00, 1.089e+00,
     &  1.086e+00, 1.086e+00, 1.087e+00, 1.088e+00, 1.088e+00,
     &  1.087e+00, 1.087e+00, 1.090e+00, 1.091e+00, 1.091e+00/
      data mrtab(651:700)/
     &  1.092e+00, 1.096e+00, 1.098e+00, 1.100e+00, 1.101e+00,
     &  1.105e+00, 1.107e+00, 1.109e+00, 1.110e+00, 1.114e+00,
     &  1.117e+00, 1.118e+00, 1.122e+00, 1.125e+00, 1.126e+00,
     &  1.131e+00, 1.134e+00, 1.136e+00, 1.140e+00, 1.141e+00,
     &  1.146e+00, 1.150e+00, 1.152e+00, 1.156e+00, 1.160e+00,
     &  1.162e+00, 1.166e+00, 1.168e+00, 1.172e+00, 1.174e+00,
     &  1.179e+00, 1.182e+00, 1.185e+00, 1.189e+00, 1.191e+00,
     &  1.196e+00, 1.198e+00, 1.202e+00, 1.204e+00, 1.207e+00,
     &  1.209e+00, 1.214e+00, 1.215e+00, 1.219e+00, 1.221e+00,
     &  1.226e+00, 1.228e+00, 1.232e+00, 1.234e+00, 1.238e+00/
      data mrtab(701:750)/
     &  1.239e+00, 1.243e+00, 1.245e+00, 1.249e+00, 1.252e+00,
     &  1.254e+00, 1.257e+00, 1.259e+00, 1.263e+00, 1.265e+00,
     &  1.268e+00, 1.270e+00, 1.273e+00, 1.276e+00, 1.278e+00,
     &  1.283e+00, 1.285e+00, 1.287e+00, 1.291e+00, 1.293e+00,
     &  1.295e+00, 1.299e+00, 1.301e+00, 1.303e+00, 1.307e+00,
     &  1.309e+00, 1.313e+00, 1.315e+00, 1.317e+00, 1.321e+00,
     &  1.323e+00, 1.325e+00, 1.327e+00, 1.331e+00, 1.333e+00,
     &  1.335e+00, 1.339e+00, 1.341e+00, 1.343e+00, 1.347e+00,
     &  1.350e+00, 1.351e+00, 1.353e+00, 1.357e+00, 1.359e+00,
     &  1.361e+00, 1.365e+00, 1.367e+00, 1.369e+00, 1.371e+00/
      data mrtab(751:800)/
     &  1.375e+00, 1.377e+00, 1.379e+00, 1.380e+00, 1.384e+00,
     &  1.386e+00, 1.388e+00, 1.389e+00, 1.391e+00, 1.394e+00,
     &  1.396e+00, 1.398e+00, 1.400e+00, 1.401e+00, 1.405e+00,
     &  1.407e+00, 1.409e+00, 1.410e+00, 1.412e+00, 1.415e+00,
     &  1.418e+00, 1.420e+00, 1.422e+00, 1.423e+00, 1.426e+00,
     &  1.428e+00, 1.430e+00, 1.431e+00, 1.433e+00, 1.435e+00,
     &  1.436e+00, 1.440e+00, 1.442e+00, 1.443e+00, 1.444e+00,
     &  1.445e+00, 1.447e+00, 1.448e+00, 1.450e+00, 1.452e+00,
     &  1.454e+00, 1.456e+00, 1.457e+00, 1.458e+00, 1.459e+00,
     &  1.460e+00, 1.460e+00, 1.463e+00, 1.465e+00, 1.467e+00/
      data mrtab(801:860)/
     &  1.467e+00, 1.467e+00, 1.468e+00, 1.484e+00, 1.499e+00,
     &  1.516e+00, 1.529e+00, 1.538e+00, 1.544e+00, 1.547e+00,
     &  1.546e+00, 1.543e+00, 1.536e+00, 1.527e+00, 1.519e+00,
     &  1.512e+00, 1.506e+00, 1.499e+00, 1.493e+00, 1.487e+00,
     &  1.481e+00, 1.478e+00, 1.477e+00, 1.476e+00, 1.477e+00,
     &  1.481e+00, 1.485e+00, 1.492e+00, 1.499e+00, 1.509e+00,
     &  1.520e+00, 1.531e+00, 1.542e+00, 1.567e+00, 1.594e+00,
     &  1.619e+00, 1.644e+00, 1.669e+00, 1.690e+00, 1.710e+00,
     &  1.729e+00, 1.747e+00, 1.763e+00, 1.777e+00, 1.791e+00,
     &  1.806e+00, 1.819e+00, 1.831e+00, 1.842e+00, 1.852e+00,
     &  1.860e+00, 1.867e+00, 1.874e+00, 1.881e+00, 1.886e+00,
     &  1.891e+00, 1.895e+00, 1.899e+00, 1.908e+00, 1.912e+00/

      data mitab(1:50)/
     &  1.745e-03, 2.370e-03, 3.146e-03, 4.072e-03, 5.174e-03,
     &  7.958e-03, 1.164e-02, 1.636e-02, 2.227e-02, 2.950e-02,
     &  3.818e-02, 4.850e-02, 6.064e-02, 7.461e-02, 9.074e-02,
     &  1.093e-01, 1.303e-01, 1.534e-01, 1.798e-01, 2.088e-01,
     &  2.414e-01, 2.766e-01, 2.998e-01, 3.154e-01, 3.310e-01,
     &  3.498e-01, 3.739e-01, 4.119e-01, 4.558e-01, 5.033e-01,
     &  5.355e-01, 5.634e-01, 5.792e-01, 5.859e-01, 5.805e-01,
     &  5.859e-01, 5.981e-01, 6.135e-01, 6.292e-01, 6.453e-01,
     &  6.573e-01, 6.573e-01, 6.528e-01, 6.439e-01, 6.292e-01,
     &  6.050e-01, 5.752e-01, 5.430e-01, 5.185e-01, 4.929e-01/
      data mitab(51:100)/
     &  4.708e-01, 4.485e-01, 4.303e-01, 4.148e-01, 3.988e-01,
     &  3.826e-01, 3.705e-01, 3.596e-01, 3.490e-01, 3.387e-01,
     &  3.220e-01, 2.876e-01, 2.392e-01, 1.870e-01, 1.489e-01,
     &  1.333e-01, 1.422e-01, 1.678e-01, 1.927e-01, 2.167e-01,
     &  2.409e-01, 2.641e-01, 2.772e-01, 2.772e-01, 2.581e-01,
     &  2.338e-01, 2.022e-01, 1.670e-01, 1.351e-01, 1.039e-01,
     &  7.241e-02, 3.997e-02, 3.998e-03, 2.004e-03, 1.182e-03,
     &  8.391e-04, 5.995e-05, 1.250e-06, 3.622e-07, 1.850e-07,
     &  1.101e-07, 6.711e-08, 3.844e-08, 1.999e-08, 1.270e-08,
     &  1.158e-08, 1.101e-08, 1.071e-08, 1.049e-08, 9.903e-09/
      data mitab(101:150)/
     &  9.307e-09, 8.606e-09, 7.994e-09, 7.444e-09, 6.852e-09,
     &  6.292e-09, 5.792e-09, 5.405e-09, 4.795e-09, 4.403e-09,
     &  4.148e-09, 3.826e-09, 3.546e-09, 3.325e-09, 3.190e-09,
     &  3.082e-09, 2.984e-09, 2.883e-09, 2.766e-09, 2.653e-09,
     &  2.528e-09, 2.420e-09, 2.316e-09, 2.217e-09, 2.117e-09,
     &  2.031e-09, 1.940e-09, 1.840e-09, 1.761e-09, 1.663e-09,
     &  1.581e-09, 1.489e-09, 1.422e-09, 1.339e-09, 1.258e-09,
     &  1.169e-09, 1.088e-09, 1.018e-09, 9.393e-10, 8.685e-10,
     &  8.087e-10, 7.795e-10, 7.600e-10, 7.495e-10, 7.291e-10,
     &  7.011e-10, 7.092e-10, 7.158e-10, 7.342e-10, 7.849e-10/
      data mitab(151:200)/
     &  9.242e-10, 1.078e-09, 1.267e-09, 1.461e-09, 1.570e-09,
     &  1.640e-09, 1.757e-09, 1.887e-09, 2.098e-09, 2.269e-09,
     &  2.442e-09, 2.659e-09, 2.869e-09, 3.132e-09, 3.434e-09,
     &  3.844e-09, 4.434e-09, 5.221e-09, 6.365e-09, 7.723e-09,
     &  9.634e-09, 1.132e-08, 1.238e-08, 1.330e-08, 1.399e-08,
     &  1.472e-08, 1.502e-08, 1.552e-08, 1.570e-08, 1.606e-08,
     &  1.674e-08, 1.777e-08, 1.940e-08, 2.031e-08, 2.098e-08,
     &  2.177e-08, 2.300e-08, 2.470e-08, 2.653e-08, 2.963e-08,
     &  3.348e-08, 4.100e-08, 4.998e-08, 5.995e-08, 7.291e-08,
     &  9.137e-08, 1.150e-07, 1.348e-07, 1.458e-07, 1.530e-07/
      data mitab(201:250)/
     &  1.559e-07, 1.580e-07, 1.580e-07, 1.570e-07, 1.527e-07,
     &  1.478e-07, 1.409e-07, 1.339e-07, 1.282e-07, 1.258e-07,
     &  1.250e-07, 1.270e-07, 1.330e-07, 1.448e-07, 1.621e-07,
     &  1.819e-07, 2.041e-07, 2.243e-07, 2.459e-07, 2.690e-07,
     &  2.929e-07, 3.153e-07, 3.348e-07, 3.546e-07, 3.748e-07,
     &  3.907e-07, 4.053e-07, 4.234e-07, 4.403e-07, 4.622e-07,
     &  4.862e-07, 5.150e-07, 5.699e-07, 6.696e-07, 8.304e-07,
     &  1.060e-06, 1.368e-06, 1.771e-06, 2.169e-06, 2.557e-06,
     &  2.932e-06, 3.190e-06, 3.358e-06, 3.464e-06, 3.502e-06,
     &  3.480e-06, 3.418e-06, 3.337e-06, 3.253e-06, 3.131e-06/
      data mitab(251:300)/
     &  3.000e-06, 2.688e-06, 2.352e-06, 2.001e-06, 1.690e-06,
     &  1.419e-06, 1.299e-06, 1.259e-06, 1.329e-06, 1.499e-06,
     &  1.708e-06, 2.038e-06, 2.628e-06, 3.869e-06, 5.951e-06,
     &  9.306e-06, 1.069e-05, 1.120e-05, 1.160e-05, 1.181e-05,
     &  1.199e-05, 1.191e-05, 1.179e-05, 1.160e-05, 1.139e-05,
     &  1.100e-05, 1.079e-05, 1.090e-05, 1.139e-05, 1.221e-05,
     &  1.400e-05, 1.639e-05, 1.912e-05, 2.251e-05, 2.849e-05,
     &  4.047e-05, 4.505e-05, 5.804e-05, 7.802e-05, 1.060e-04,
     &  1.530e-04, 2.540e-04, 3.197e-04, 3.538e-04, 3.629e-04,
     &  3.637e-04, 3.604e-04, 3.387e-04, 3.018e-04, 2.659e-04/
      data mitab(301:350)/
     &  2.248e-04, 1.958e-04, 1.741e-04, 1.602e-04, 1.441e-04,
     &  1.348e-04, 1.240e-04, 1.140e-04, 1.071e-04, 9.940e-05,
     &  9.347e-05, 8.804e-05, 8.310e-05, 8.096e-05, 7.903e-05,
     &  7.591e-05, 7.398e-05, 7.404e-05, 7.495e-05, 7.601e-05,
     &  7.743e-05, 8.050e-05, 8.410e-05, 8.900e-05, 9.510e-05,
     &  1.000e-04, 1.051e-04, 1.120e-04, 1.219e-04, 1.330e-04,
     &  1.359e-04, 1.371e-04, 1.380e-04, 1.418e-04, 1.552e-04,
     &  1.861e-04, 3.205e-04, 5.209e-04, 7.224e-04, 9.221e-04,
     &  1.161e-03, 1.678e-03, 1.827e-03, 1.922e-03, 1.909e-03,
     &  1.848e-03, 1.717e-03, 1.548e-03, 1.402e-03, 1.250e-03/
      data mitab(351:400)/
     &  1.101e-03, 9.904e-04, 8.888e-04, 8.050e-04, 7.392e-04,
     &  6.742e-04, 6.206e-04, 5.725e-04, 5.294e-04, 4.884e-04,
     &  4.643e-04, 4.403e-04, 4.176e-04, 3.970e-04, 3.826e-04,
     &  3.705e-04, 3.588e-04, 3.506e-04, 3.434e-04, 3.395e-04,
     &  3.379e-04, 3.387e-04, 3.410e-04, 3.458e-04, 3.571e-04,
     &  3.739e-04, 3.898e-04, 4.081e-04, 4.293e-04, 4.506e-04,
     &  4.686e-04, 4.918e-04, 5.114e-04, 5.430e-04, 5.995e-04,
     &  6.365e-04, 6.852e-04, 7.427e-04, 7.921e-04, 8.488e-04,
     &  9.095e-04, 9.904e-04, 1.071e-03, 1.150e-03, 1.250e-03,
     &  1.348e-03, 1.472e-03, 1.581e-03, 1.709e-03, 1.811e-03/
      data mitab(401:450)/
     &  1.900e-03, 1.953e-03, 1.990e-03, 2.017e-03, 2.069e-03,
     &  2.142e-03, 2.269e-03, 2.311e-03, 2.338e-03, 2.387e-03,
     &  2.425e-03, 2.476e-03, 2.575e-03, 2.703e-03, 2.977e-03,
     &  3.302e-03, 4.016e-03, 4.363e-03, 4.828e-03, 5.368e-03,
     &  6.278e-03, 7.325e-03, 8.547e-03, 1.049e-02, 1.270e-02,
     &  1.451e-02, 1.640e-02, 1.861e-02, 2.050e-02, 2.817e-02,
     &  3.800e-02, 4.622e-02, 5.480e-02, 6.483e-02, 7.444e-02,
     &  8.352e-02, 9.285e-02, 1.020e-01, 1.119e-01, 1.210e-01,
     &  1.312e-01, 1.422e-01, 1.541e-01, 1.670e-01, 1.798e-01,
     &  1.940e-01, 2.060e-01, 2.182e-01, 2.290e-01, 2.392e-01/
      data mitab(451:500)/
     &  2.493e-01, 2.581e-01, 2.647e-01, 2.715e-01, 2.759e-01,
     &  2.798e-01, 2.804e-01, 2.823e-01, 2.817e-01, 2.785e-01,
     &  2.759e-01, 2.721e-01, 2.721e-01, 2.398e-01, 1.918e-01,
     &  1.348e-01, 9.242e-02, 6.107e-02, 3.688e-02, 2.611e-02,
     &  1.949e-02, 1.321e-02, 9.393e-03, 6.789e-03, 5.150e-03,
     &  4.234e-03, 3.596e-03, 3.402e-03, 3.402e-03, 3.530e-03,
     &  3.800e-03, 4.157e-03, 4.600e-03, 5.067e-03, 5.621e-03,
     &  6.220e-03, 6.883e-03, 7.600e-03, 8.449e-03, 9.307e-03,
     &  1.030e-02, 1.140e-02, 1.238e-02, 1.361e-02, 1.472e-02,
     &  1.548e-02, 1.570e-02, 1.552e-02, 1.499e-02, 1.441e-02/
      data mitab(501:550)/
     &  1.370e-02, 1.312e-02, 1.241e-02, 1.180e-02, 1.111e-02,
     &  1.061e-02, 1.011e-02, 9.904e-03, 9.790e-03, 9.881e-03,
     &  1.030e-02, 1.078e-02, 1.158e-02, 1.258e-02, 1.418e-02,
     &  1.659e-02, 2.031e-02, 2.482e-02, 3.295e-02, 4.323e-02,
     &  6.220e-02, 8.646e-02, 1.069e-01, 1.250e-01, 1.308e-01,
     &  1.172e-01, 8.786e-02, 6.947e-02, 5.699e-02, 4.952e-02,
     &  4.485e-02, 4.176e-02, 3.925e-02, 3.731e-02, 3.563e-02,
     &  3.450e-02, 3.371e-02, 3.310e-02, 3.272e-02, 3.242e-02,
     &  3.220e-02, 3.212e-02, 3.197e-02, 3.190e-02, 3.197e-02,
     &  3.205e-02, 3.205e-02, 3.220e-02, 3.220e-02, 3.227e-02/
      data mitab(551:600)/
     &  3.242e-02, 3.249e-02, 3.257e-02, 3.272e-02, 3.279e-02,
     &  3.302e-02, 3.310e-02, 3.325e-02, 3.348e-02, 3.371e-02,
     &  3.395e-02, 3.410e-02, 3.426e-02, 3.450e-02, 3.466e-02,
     &  3.490e-02, 3.514e-02, 3.530e-02, 3.563e-02, 3.579e-02,
     &  3.604e-02, 3.637e-02, 3.654e-02, 3.688e-02, 3.714e-02,
     &  3.748e-02, 3.783e-02, 3.809e-02, 3.844e-02, 3.880e-02,
     &  3.916e-02, 3.952e-02, 3.988e-02, 4.034e-02, 4.072e-02,
     &  4.109e-02, 4.147e-02, 4.196e-02, 4.234e-02, 4.293e-02,
     &  4.333e-02, 4.373e-02, 4.434e-02, 4.475e-02, 4.537e-02,
     &  4.600e-02, 4.664e-02, 4.718e-02, 4.784e-02, 4.850e-02/
      data mitab(601:650)/
     &  4.929e-02, 4.998e-02, 5.079e-02, 5.174e-02, 5.270e-02,
     &  5.380e-02, 5.805e-02, 5.634e-02, 5.845e-02, 5.995e-02,
     &  6.191e-02, 6.394e-02, 6.619e-02, 6.852e-02, 7.092e-02,
     &  7.358e-02, 7.652e-02, 7.958e-02, 8.295e-02, 8.646e-02,
     &  8.970e-02, 9.328e-02, 9.678e-02, 9.995e-02, 1.039e-01,
     &  1.083e-01, 1.129e-01, 1.172e-01, 1.218e-01, 1.270e-01,
     &  1.321e-01, 1.370e-01, 1.422e-01, 1.472e-01, 1.520e-01,
     &  1.570e-01, 1.621e-01, 1.678e-01, 1.741e-01, 1.802e-01,
     &  1.865e-01, 1.927e-01, 1.990e-01, 2.055e-01, 2.112e-01,
     &  2.177e-01, 2.238e-01, 2.295e-01, 2.359e-01, 2.420e-01/
      data mitab(651:700)/
     &  2.476e-01, 2.528e-01, 2.593e-01, 2.641e-01, 2.690e-01,
     &  2.740e-01, 2.791e-01, 2.837e-01, 2.883e-01, 2.929e-01,
     &  2.977e-01, 3.012e-01, 3.060e-01, 3.103e-01, 3.139e-01,
     &  3.183e-01, 3.227e-01, 3.257e-01, 3.295e-01, 3.325e-01,
     &  3.363e-01, 3.402e-01, 3.426e-01, 3.466e-01, 3.490e-01,
     &  3.514e-01, 3.546e-01, 3.571e-01, 3.596e-01, 3.621e-01,
     &  3.646e-01, 3.679e-01, 3.696e-01, 3.722e-01, 3.739e-01,
     &  3.757e-01, 3.774e-01, 3.791e-01, 3.809e-01, 3.826e-01,
     &  3.844e-01, 3.862e-01, 3.871e-01, 3.897e-01, 3.916e-01,
     &  3.934e-01, 3.943e-01, 3.961e-01, 3.970e-01, 3.988e-01/
      data mitab(701:750)/
     &  3.997e-01, 4.016e-01, 4.025e-01, 4.044e-01, 4.053e-01,
     &  4.053e-01, 4.072e-01, 4.081e-01, 4.100e-01, 4.100e-01,
     &  4.119e-01, 4.128e-01, 4.138e-01, 4.148e-01, 4.157e-01,
     &  4.176e-01, 4.176e-01, 4.186e-01, 4.196e-01, 4.205e-01,
     &  4.205e-01, 4.225e-01, 4.225e-01, 4.225e-01, 4.234e-01,
     &  4.244e-01, 4.254e-01, 4.254e-01, 4.254e-01, 4.264e-01,
     &  4.264e-01, 4.274e-01, 4.274e-01, 4.283e-01, 4.283e-01,
     &  4.283e-01, 4.293e-01, 4.293e-01, 4.293e-01, 4.303e-01,
     &  4.293e-01, 4.293e-01, 4.293e-01, 4.303e-01, 4.303e-01,
     &  4.293e-01, 4.303e-01, 4.303e-01, 4.293e-01, 4.293e-01/
      data mitab(751:800)/
     &  4.293e-01, 4.293e-01, 4.283e-01, 4.283e-01, 4.283e-01,
     &  4.283e-01, 4.274e-01, 4.274e-01, 4.264e-01, 4.274e-01,
     &  4.264e-01, 4.264e-01, 4.254e-01, 4.254e-01, 4.254e-01,
     &  4.254e-01, 4.244e-01, 4.244e-01, 4.234e-01, 4.234e-01,
     &  4.234e-01, 4.225e-01, 4.215e-01, 4.205e-01, 4.205e-01,
     &  4.196e-01, 4.186e-01, 4.176e-01, 4.176e-01, 4.167e-01,
     &  4.157e-01, 4.157e-01, 4.138e-01, 4.128e-01, 4.119e-01,
     &  4.109e-01, 4.100e-01, 4.091e-01, 4.091e-01, 4.081e-01,
     &  4.072e-01, 4.062e-01, 4.044e-01, 4.034e-01, 4.025e-01,
     &  4.016e-01, 4.007e-01, 4.007e-01, 3.988e-01, 3.970e-01/
      data mitab(801:860)/
     &  3.952e-01, 3.943e-01, 3.934e-01, 3.818e-01, 3.722e-01,
     &  3.629e-01, 3.482e-01, 3.356e-01, 3.227e-01, 3.103e-01,
     &  2.991e-01, 2.889e-01, 2.817e-01, 2.791e-01, 2.798e-01,
     &  2.830e-01, 2.863e-01, 2.916e-01, 2.991e-01, 3.068e-01,
     &  3.190e-01, 3.317e-01, 3.442e-01, 3.588e-01, 3.739e-01,
     &  3.880e-01, 4.025e-01, 4.176e-01, 4.323e-01, 4.465e-01,
     &  4.579e-01, 4.675e-01, 4.773e-01, 4.975e-01, 5.079e-01,
     &  5.162e-01, 5.233e-01, 5.258e-01, 5.245e-01, 5.245e-01,
     &  5.233e-01, 5.209e-01, 5.174e-01, 5.138e-01, 5.103e-01,
     &  5.079e-01, 5.021e-01, 4.964e-01, 4.907e-01, 4.839e-01,
     &  4.773e-01, 4.718e-01, 4.654e-01, 4.600e-01, 4.548e-01,
     &  4.485e-01, 4.434e-01, 4.383e-01, 4.176e-01, 4.167e-01/

      call locate(wltab,ntab,wl,i)
      wt=(wl-wltab(i))/(wltab(i+1)-wltab(i))
      wt=max(0.0_kr,min(1.0_kr,wt))     ! no extrapolation
      nr=mrtab(i)*(mrtab(i+1)/mrtab(i))**wt
      ni=mitab(i)*(mitab(i+1)/mitab(i))**wt

c 
c correction to be applied to the index of refraction and to the extinction 
c coefficients of the pure water to obtain the ocean water one (see for 
c example friedman). by default, a typical sea water is assumed 
c (salinity=34.3ppt, chlorinity=19ppt) as reported by sverdrup. 
c in that case there is no correction for the extinction coefficient between 
c 0.25 and 4 microns. for the index of refraction, a correction of +0.006 
c has to be applied (mclellan). for a chlorinity of 19.0ppt the correction 
c is a linear function of the salt concentration. 

c references:
c friedman d., applied optics, 1969, vol.8, no.10, pp.2073-2078.
c mclellan h.j., elements of physical oceanography, pergamon press, inc.,
c        new-york, 1965, p 129.
c sverdrup h.v. et al., the oceans (prentice-hall, inc., englewood cliffs,
c        n.j., 1942, p 173.

        nr=nr+nrc*(xsal/34.3)
        ni=ni+nic*(xsal/34.3)

        return
        end
c
c=======================================================================
      subroutine sunglint(wndspd,nr,ni,csin,cvin,phi,rog)

c input parameters:   wndspd=speed of the wind (in m/s)
c                     nr=index of refraction of the sea water
c                     ni=extinction coefficient of the sea water
c                     cs=cosine of solar zenith angle
c                     cv=cosine of view zenith angle
c                     phi=relative azimuth (sun-satellite) (in radians)
c
c output parameters:  rog=reflectance of the sun glint
c
      use params, only: kr, pi, zip
      implicit none

      real(kr), intent(in) :: wndspd,nr,ni,csin,cvin,phi
      real(kr), intent(out):: rog

      real(kr) ::  d2r,phw,cs,cv,ss,sv,zx,zy,tantilt,tilt,proba,
     &     xe,xn,xe2,xn2,coef,cos2chi,coschi,sinchi,r1,sigmac,
     &     sigmau,c21,c03,c40,c04,c22,axe2,axn2,axe4,axn4,axe2xn2

      d2r=pi/180.

c  no wave shodowing, so set a floor on incident and viewing angle 

      cs=max(csin,0.05_kr) 
      cv=max(cvin,0.05_kr) 
      ss=sqrt(1.-cs**2)
      sv=sqrt(1.-cv**2)

      zx=-sv*sin(pi-phi)/(cs+cv)
      zy=(ss+sv*cos(pi-phi))/(cs+cv)
      tantilt=sqrt(zx*zx+zy*zy)
      tilt=atan(tantilt)

c  anisotropic gaussian distribution

      sigmac=0.003+0.00192*wndspd
      sigmau=0.00316*wndspd
      c21=0.01-0.0086*wndspd
      c03=0.04-0.033*wndspd
      c40=0.40
      c22=0.12
      c04=0.23

                                ! brdf dependence on wndang, 
                                ! original code from 6s:
                                
c     phw=wndang*d2r            ! phw=phi_sun-phi_wind
c     xe=(cos(phw)*zx+sin(phw)*zy)/sqrt(sigmac)
c     xn=(-sin(phw)*zx+cos(phw)*zy)/sqrt(sigmau)
c     xe2=xe*xe
c     xn2=xn*xn
c     coef=1-c21/2.*(xe2-1)*xn-c03/6.*(xn2-3)*xn
c     coef=coef+c40/24.*(xe2*xe2-6*xe2+3)
c     coef=coef+c04/24.*(xn2*xn2-6*xn2+3)
c     coef=coef+c22/4.*(xe2-1)*(xn2-1)
c     coef=coef/(2.*pi*sqrt(sigmau)*sqrt(sigmac))
c     proba=coef*exp(-(xe2+xn2)/2.)

                                ! angle average over wind directions
                                ! axe2=<xe**2> averaged over all angles phw
                                ! <xe**3>=<xn**3>=0
                                ! <xe>=<xn>=<xn2*xe>=<xe2*xn>=0
                                ! <sin2>=1/2 <sin4>=3/8 <sin2 cos2>=1/8

      axe2=.5*(zx**2+zy**2)/sigmac
      axn2=.5*(zx**2+zy**2)/sigmau
      axe4=(3*zx**4+6*zx**2*zy**2+3*zy**4)/(8*sigmac**2)
      axn4=(3*zx**4+6*zx**2*zy**2+3*zy**4)/(8*sigmau**2)
      axe2xn2=(zx**4+10*zx**2*zy**2+zy**4)/(8*sigmau*sigmac)        
      coef=1.
      coef=coef+c40/24.*(axe4-6*axe2+3)
      coef=coef+c04/24.*(axn4-6*axn2+3)
      coef=coef+c22/4.*(axe2xn2-axn2-axe2+1)
      coef=coef/(2.*pi*sqrt(sigmau)*sqrt(sigmac))
      proba=coef*exp(-(axe2+axn2)/2.)

c compute fresnel's coefficient r1

      cos2chi=cv*cs+sv*ss*cos(pi-phi)
      if (cos2chi.gt.1.0)cos2chi=0.99999999999
      if (cos2chi.lt.-1.0)cos2chi=-0.99999999999
      coschi=sqrt(0.5*(1+cos2chi))
      sinchi=sqrt(0.5*(1-cos2chi))
      call fresnel(nr,ni,coschi,sinchi,r1)

c compute reflectance of the sun glint

      rog=pi*r1*proba/(4.*cs*cv*cos(tilt)**4)

      return
      end
c
c
c=======================================================================
      subroutine fresnel(nr,ni,coschi,sinchi,r1)
c
c to compute the fresnel's coefficient of reflection (see for
c example m. born and e. wolf, principles of optics, pergamon press, fifth
c edition, 1975, pp 628
c input parameters: nr=index of refraction of the sea water
c                   ni=extinction coefficient of the sea water
c                   coschi & sinchi=cosine and sine of the incident radiation 
c                                   with respect of the wave facet normal.
c output parameter: r1=fresnel's coefficient for reflection (intensity)
c
c rr2=intensity reflection coefficient for perpendicular component
c rl2=intensity reflection coefficient for parallel component
c r1=total unpolarized reflection coefficient

      use params, only: kr
      implicit none
      real(kr) ::  nr,ni,a1,a2,u,v,rr2,rl2,b1,b2,r1,coschi,sinchi

c absolute value for a1 to get v=0 when ni=0

      a1=abs(nr*nr-ni*ni-sinchi*sinchi)
      a2=sqrt((nr*nr-ni*ni-sinchi*sinchi)**2.+4*nr*nr*ni*ni)
      u=sqrt(0.5*(a1+a2))
      v=sqrt(0.5*(-a1+a2))
      rr2=((coschi-u)**2+v*v)/((coschi+u)**2+v*v)
      b1=(nr*nr-ni*ni)*coschi
      b2=2*nr*ni*coschi
      rl2=((b1-u)**2+(b2+v)**2)/((b1+u)**2+(b2-v)**2)

      r1=(rr2+rl2)/2.

      return
      end
c
c
c=======================================================================

      module sunblk             ! solar data base
      use params, only: kr
      implicit none
      save
      integer, parameter :: mxwl=5000
      real(kr), dimension(mxwl) :: wlsun,sun
      integer :: nns=0
      end module sunblk
c=======================================================================
      function solirr(wl,nf)
c
c input:
c     wl         wavelength in microns
c
c     nf         solar spectrum database switch
c               -1 = read solar.dat, replace data in wlsun sun
c                0 = spectrally uniform
c                1 = 5s 
c                2 = lowtran7 
c                3 = modtran3 (degraded to 20 cm-1 resolution)
c
c    note: the value of nf is not used after the first call to solirr
c    note: the 5s model uses a 1/w**4 power law for w > 4 um
c
c output:
c                extra-terestrial solar flux in w/m2/micron
c
      use params, only: zero, one, kr
      use sunblk, only: mxwl, nns, wlsun, sun
      implicit none
      integer :: j, nf
      real(kr) :: wt, wl, solirr

      if(nf.eq.0) then
        solirr=1.
        return
      endif

      if(nns.eq.0) then
        nns=mxwl
        select case (nf)
        case (-1) ; call rdspec('solar.dat',wlsun,sun,nns)
        case (1)  ; call sun1s(wlsun,sun,nns)
        case (2)  ; call sunlow(wlsun,sun,nns)
        case (3)  ; call sunmod(wlsun,sun,nns)
        case default 
           print *,'ERROR in solirr --- illegal value of nf',nf
           stop
        end select
      endif

      call locate(wlsun,nns,wl,j)
      wt=(wl-wlsun(j))/(wlsun(j+1)-wlsun(j))
      wt=max(zero,min(one,wt))
      solirr=sun(j)*(1.-wt)+sun(j+1)*wt

      return
      end      
c=======================================================================
      subroutine sun1s(wlsun,sun,nns)
      use params, only: kr
      implicit none
      integer :: i, nns
      real(kr) :: wlsun(*),sun(*),sun1(751)
c
c block    source           range           resolution      n_elements
c -----    ------           -----           ----------      ----------
c sun1     5s           0.25-4.00 um        .005 um            751
c
c units are w/m2/micron
c

      data sun1(1:60)/
     &  6.930e+01,8.600e+01,1.141e+02,1.600e+02,1.790e+02,1.770e+02,
     &  2.104e+02,2.730e+02,3.638e+02,5.050e+02,5.580e+02,5.370e+02,
     &  5.810e+02,6.580e+02,7.308e+02,8.180e+02,9.260e+02,8.981e+02,
     &  9.433e+02,9.087e+02,9.635e+02,1.026e+03,9.351e+02,1.137e+03,
     &  1.190e+03,1.029e+03,1.149e+03,9.181e+02,1.146e+03,9.245e+02,
     &  1.614e+03,1.648e+03,1.706e+03,1.783e+03,1.717e+03,1.693e+03,
     &  1.502e+03,1.760e+03,1.790e+03,1.929e+03,2.057e+03,2.021e+03,
     &  2.040e+03,2.013e+03,1.986e+03,2.019e+03,2.057e+03,1.879e+03,
     &  1.933e+03,1.974e+03,1.897e+03,1.936e+03,1.938e+03,1.825e+03,
     &  1.803e+03,1.861e+03,1.911e+03,1.899e+03,1.853e+03,1.878e+03/
      data sun1(61:120)/
     &  1.871e+03,1.874e+03,1.827e+03,1.847e+03,1.835e+03,1.864e+03,
     &  1.838e+03,1.843e+03,1.763e+03,1.795e+03,1.760e+03,1.769e+03,
     &  1.737e+03,1.689e+03,1.717e+03,1.675e+03,1.669e+03,1.658e+03,
     &  1.637e+03,1.623e+03,1.589e+03,1.515e+03,1.543e+03,1.555e+03,
     &  1.534e+03,1.517e+03,1.499e+03,1.469e+03,1.459e+03,1.438e+03,
     &  1.408e+03,1.407e+03,1.389e+03,1.368e+03,1.336e+03,1.352e+03,
     &  1.328e+03,1.314e+03,1.276e+03,1.286e+03,1.270e+03,1.257e+03,
     &  1.242e+03,1.220e+03,1.204e+03,1.199e+03,1.190e+03,1.181e+03,
     &  1.168e+03,1.144e+03,1.143e+03,1.123e+03,1.110e+03,1.112e+03,
     &  1.078e+03,1.077e+03,1.070e+03,1.047e+03,1.044e+03,1.030e+03/
      data sun1(121:180)/
     &  9.909e+02,9.429e+02,1.003e+03,9.550e+02,9.669e+02,9.997e+02,
     &  9.841e+02,9.740e+02,9.636e+02,9.529e+02,9.430e+02,9.330e+02,
     &  9.230e+02,9.134e+02,9.041e+02,8.940e+02,8.844e+02,8.751e+02,
     &  8.654e+02,8.565e+02,8.475e+02,8.385e+02,8.295e+02,8.205e+02,
     &  8.115e+02,8.025e+02,7.939e+02,7.857e+02,7.769e+02,7.690e+02,
     &  7.606e+02,7.515e+02,7.425e+02,7.339e+02,7.256e+02,7.169e+02,
     &  7.090e+02,7.006e+02,6.919e+02,6.840e+02,6.760e+02,6.676e+02,
     &  6.589e+02,6.518e+02,6.464e+02,6.411e+02,6.350e+02,6.294e+02,
     &  6.241e+02,6.180e+02,6.124e+02,6.075e+02,6.021e+02,5.964e+02,
     &  5.915e+02,5.861e+02,5.804e+02,5.755e+02,5.705e+02,5.655e+02/
      data sun1(181:240)/
     &  5.605e+02,5.505e+02,5.505e+02,5.455e+02,5.405e+02,5.355e+02,
     &  5.309e+02,5.270e+02,5.226e+02,5.179e+02,5.130e+02,5.062e+02,
     &  4.996e+02,4.931e+02,4.865e+02,4.800e+02,4.735e+02,4.669e+02,
     &  4.604e+02,4.538e+02,4.473e+02,4.408e+02,4.342e+02,4.277e+02,
     &  4.211e+02,4.146e+02,4.080e+02,4.015e+02,3.950e+02,3.884e+02,
     &  3.967e+02,4.071e+02,4.035e+02,3.998e+02,3.962e+02,3.925e+02,
     &  3.889e+02,3.852e+02,3.816e+02,3.780e+02,3.743e+02,3.707e+02,
     &  3.670e+02,3.634e+02,3.597e+02,3.560e+02,3.524e+02,3.488e+02,
     &  3.451e+02,3.415e+02,3.422e+02,3.437e+02,3.409e+02,3.381e+02,
     &  3.353e+02,3.325e+02,3.297e+02,3.270e+02,3.242e+02,3.214e+02/
      data sun1(241:300)/
     &  3.186e+02,3.158e+02,3.130e+02,3.102e+02,3.075e+02,3.047e+02,
     &  3.019e+02,2.991e+02,2.963e+02,2.935e+02,2.933e+02,2.934e+02,
     &  2.911e+02,2.888e+02,2.865e+02,2.843e+02,2.820e+02,2.797e+02,
     &  2.774e+02,2.751e+02,2.728e+02,2.705e+02,2.682e+02,2.659e+02,
     &  2.636e+02,2.613e+02,2.591e+02,2.568e+02,2.545e+02,2.522e+02,
     &  2.501e+02,2.481e+02,2.459e+02,2.436e+02,2.414e+02,2.391e+02,
     &  2.369e+02,2.346e+02,2.324e+02,2.301e+02,2.279e+02,2.256e+02,
     &  2.234e+02,2.212e+02,2.189e+02,2.167e+02,2.144e+02,2.122e+02,
     &  2.099e+02,2.077e+02,2.063e+02,2.050e+02,2.029e+02,2.009e+02,
     &  1.988e+02,1.967e+02,1.946e+02,1.925e+02,1.905e+02,1.884e+02/
      data sun1(301:360)/
     &  1.863e+02,1.842e+02,1.821e+02,1.801e+02,1.780e+02,1.759e+02,
     &  1.738e+02,1.717e+02,1.697e+02,1.676e+02,1.676e+02,1.679e+02,
     &  1.663e+02,1.646e+02,1.629e+02,1.613e+02,1.596e+02,1.579e+02,
     &  1.562e+02,1.546e+02,1.529e+02,1.512e+02,1.496e+02,1.479e+02,
     &  1.462e+02,1.446e+02,1.429e+02,1.412e+02,1.395e+02,1.379e+02,
     &  1.381e+02,1.386e+02,1.373e+02,1.360e+02,1.347e+02,1.334e+02,
     &  1.321e+02,1.308e+02,1.295e+02,1.282e+02,1.269e+02,1.256e+02,
     &  1.243e+02,1.230e+02,1.217e+02,1.204e+02,1.191e+02,1.178e+02,
     &  1.165e+02,1.152e+02,1.153e+02,1.156e+02,1.146e+02,1.136e+02,
     &  1.125e+02,1.115e+02,1.105e+02,1.095e+02,1.084e+02,1.074e+02/
      data sun1(361:420)/
     &  1.064e+02,1.054e+02,1.043e+02,1.033e+02,1.023e+02,1.013e+02,
     &  1.003e+02,9.922e+01,9.820e+01,9.718e+01,9.707e+01,9.711e+01,
     &  9.626e+01,9.542e+01,9.456e+01,9.372e+01,9.289e+01,9.203e+01,
     &  9.119e+01,9.034e+01,8.949e+01,8.866e+01,8.781e+01,8.697e+01,
     &  8.612e+01,8.526e+01,8.443e+01,8.359e+01,8.275e+01,8.189e+01,
     &  8.189e+01,8.202e+01,8.135e+01,8.065e+01,7.999e+01,7.930e+01,
     &  7.861e+01,7.794e+01,7.726e+01,7.657e+01,7.590e+01,7.522e+01,
     &  7.454e+01,7.386e+01,7.318e+01,7.250e+01,7.182e+01,7.114e+01,
     &  7.047e+01,6.979e+01,6.973e+01,6.980e+01,6.924e+01,6.868e+01,
     &  6.814e+01,6.757e+01,6.702e+01,6.648e+01,6.590e+01,6.536e+01/
      data sun1(421:480)/
     &  6.480e+01,6.425e+01,6.369e+01,6.314e+01,6.257e+01,6.203e+01,
     &  6.147e+01,6.092e+01,6.036e+01,5.981e+01,5.980e+01,5.985e+01,
     &  5.940e+01,5.895e+01,5.850e+01,5.806e+01,5.760e+01,5.715e+01,
     &  5.670e+01,5.626e+01,5.579e+01,5.536e+01,5.490e+01,5.443e+01,
     &  5.400e+01,5.355e+01,5.309e+01,5.265e+01,5.220e+01,5.175e+01,
     &  5.168e+01,5.167e+01,5.130e+01,5.092e+01,5.055e+01,5.018e+01,
     &  4.979e+01,4.943e+01,4.904e+01,4.867e+01,4.830e+01,4.793e+01,
     &  4.754e+01,4.718e+01,4.679e+01,4.642e+01,4.606e+01,4.567e+01,
     &  4.530e+01,4.493e+01,4.482e+01,4.479e+01,4.445e+01,4.414e+01,
     &  4.383e+01,4.349e+01,4.318e+01,4.286e+01,4.255e+01,4.221e+01/
      data sun1(481:540)/
     &  4.190e+01,4.158e+01,4.126e+01,4.094e+01,4.062e+01,4.031e+01,
     &  3.997e+01,3.966e+01,3.934e+01,3.903e+01,3.896e+01,3.892e+01,
     &  3.866e+01,3.838e+01,3.812e+01,3.786e+01,3.758e+01,3.730e+01,
     &  3.705e+01,3.677e+01,3.650e+01,3.623e+01,3.596e+01,3.569e+01,
     &  3.542e+01,3.515e+01,3.488e+01,3.461e+01,3.434e+01,3.407e+01,
     &  3.403e+01,3.403e+01,3.379e+01,3.357e+01,3.335e+01,3.312e+01,
     &  3.290e+01,3.267e+01,3.244e+01,3.223e+01,3.200e+01,3.177e+01,
     &  3.155e+01,3.131e+01,3.110e+01,3.087e+01,3.066e+01,3.042e+01,
     &  3.020e+01,2.997e+01,2.990e+01,2.987e+01,2.966e+01,2.946e+01,
     &  2.927e+01,2.908e+01,2.888e+01,2.867e+01,2.849e+01,2.830e+01/
      data sun1(541:600)/
     &  2.810e+01,2.791e+01,2.771e+01,2.752e+01,2.731e+01,2.712e+01,
     &  2.693e+01,2.674e+01,2.654e+01,2.635e+01,2.631e+01,2.628e+01,
     &  2.612e+01,2.595e+01,2.579e+01,2.564e+01,2.545e+01,2.530e+01,
     &  2.512e+01,2.498e+01,2.480e+01,2.463e+01,2.447e+01,2.431e+01,
     &  2.414e+01,2.397e+01,2.381e+01,2.366e+01,2.347e+01,2.331e+01,
     &  2.326e+01,2.321e+01,2.306e+01,2.292e+01,2.276e+01,2.263e+01,
     &  2.247e+01,2.235e+01,2.219e+01,2.204e+01,2.190e+01,2.176e+01,
     &  2.160e+01,2.147e+01,2.131e+01,2.118e+01,2.103e+01,2.089e+01,
     &  2.073e+01,2.060e+01,2.055e+01,2.053e+01,2.040e+01,2.028e+01,
     &  2.014e+01,2.003e+01,1.990e+01,1.977e+01,1.965e+01,1.953e+01/
      data sun1(601:660)/
     &  1.939e+01,1.927e+01,1.915e+01,1.902e+01,1.890e+01,1.878e+01,
     &  1.864e+01,1.853e+01,1.840e+01,1.827e+01,1.825e+01,1.824e+01,
     &  1.814e+01,1.803e+01,1.793e+01,1.783e+01,1.771e+01,1.762e+01,
     &  1.750e+01,1.741e+01,1.729e+01,1.721e+01,1.708e+01,1.700e+01,
     &  1.687e+01,1.679e+01,1.666e+01,1.657e+01,1.646e+01,1.636e+01,
     &  1.633e+01,1.630e+01,1.622e+01,1.613e+01,1.604e+01,1.596e+01,
     &  1.584e+01,1.578e+01,1.568e+01,1.559e+01,1.550e+01,1.540e+01,
     &  1.532e+01,1.524e+01,1.513e+01,1.505e+01,1.496e+01,1.487e+01,
     &  1.478e+01,1.469e+01,1.464e+01,1.463e+01,1.455e+01,1.445e+01,
     &  1.438e+01,1.432e+01,1.421e+01,1.415e+01,1.405e+01,1.398e+01/
      data sun1(661:720)/
     &  1.391e+01,1.382e+01,1.374e+01,1.367e+01,1.358e+01,1.350e+01,
     &  1.343e+01,1.335e+01,1.325e+01,1.318e+01,1.316e+01,1.312e+01,
     &  1.306e+01,1.300e+01,1.292e+01,1.285e+01,1.278e+01,1.270e+01,
     &  1.265e+01,1.257e+01,1.250e+01,1.243e+01,1.236e+01,1.228e+01,
     &  1.222e+01,1.215e+01,1.207e+01,1.201e+01,1.194e+01,1.186e+01,
     &  1.185e+01,1.184e+01,1.178e+01,1.172e+01,1.166e+01,1.160e+01,
     &  1.155e+01,1.147e+01,1.142e+01,1.136e+01,1.130e+01,1.124e+01,
     &  1.118e+01,1.112e+01,1.106e+01,1.100e+01,1.094e+01,1.089e+01,
     &  1.082e+01,1.075e+01,1.072e+01,1.070e+01,1.064e+01,1.059e+01,
     &  1.052e+01,1.047e+01,1.042e+01,1.037e+01,1.031e+01,1.025e+01/
      data sun1(721:751)/
     &  1.020e+01,1.015e+01,1.010e+01,1.003e+01,9.980e+00,9.920e+00,
     &  9.860e+00,9.820e+00,9.750e+00,9.710e+00,9.700e+00,9.700e+00,
     &  9.630e+00,9.590e+00,9.560e+00,9.500e+00,9.450e+00,9.410e+00,
     &  9.360e+00,9.320e+00,9.270e+00,9.220e+00,9.180e+00,9.130e+00,
     &  9.080e+00,9.050e+00,8.990e+00,8.940e+00,8.900e+00,8.850e+00,
     &  8.810e+00/

      nns=751
      sun(1:nns)=sun1(1:nns)
      do i=1,nns 
        wlsun(i)=.25+(4.0-.25)*real(i-1)/(nns-1)
      enddo
      return
      end
c=======================================================================
      subroutine sunlow(wlsun,sun,nns)
      use params, only: one,zero,kr
      implicit none
      integer :: i, nn1, nn2, nns
      real(kr) :: wn, wnmax, wnmin,
     &     wlsun(*), sun(*), sun2a(1440), sun2b(2910)

c block    source           range           resolution      n_elements
c -----    ------           -----           ----------      ----------
c sun2a   lowtran7        0-28780 cm-1      20 cm-1           1440
c sun2b   lowtran7    28400-57490 cm-1      10 cm-1           2910

c
      data sun2a(1:60)/
     &  0.000e+00,4.576e-08,7.010e-07,3.458e-06,1.073e-05,2.570e-05,
     &  5.250e-05,9.600e-05,1.619e-04,2.577e-04,3.910e-04,5.692e-04,
     &  8.020e-04,1.101e-03,1.477e-03,1.946e-03,2.521e-03,3.215e-03,
     &  4.044e-03,5.023e-03,6.170e-03,7.515e-03,9.068e-03,1.085e-02,
     &  1.289e-02,1.521e-02,1.776e-02,2.064e-02,2.389e-02,2.752e-02,
     &  3.154e-02,3.596e-02,4.085e-02,4.624e-02,5.213e-02,5.854e-02,
     &  6.549e-02,7.302e-02,8.117e-02,9.000e-02,9.954e-02,1.098e-01,
     &  1.208e-01,1.326e-01,1.452e-01,1.586e-01,1.731e-01,1.885e-01,
     &  2.049e-01,2.224e-01,2.411e-01,2.609e-01,2.820e-01,3.043e-01,
     &  3.279e-01,3.527e-01,3.789e-01,4.065e-01,4.355e-01,4.660e-01/
      data sun2a(61:120)/
     &  4.980e-01,5.316e-01,5.669e-01,6.039e-01,6.426e-01,6.832e-01,
     &  7.256e-01,7.699e-01,8.162e-01,8.644e-01,9.147e-01,9.671e-01,
     &  1.022e+00,1.078e+00,1.137e+00,1.199e+00,1.263e+00,1.329e+00,
     &  1.399e+00,1.471e+00,1.546e+00,1.625e+00,1.706e+00,1.791e+00,
     &  1.880e+00,1.971e+00,2.067e+00,2.166e+00,2.268e+00,2.374e+00,
     &  2.484e+00,2.597e+00,2.714e+00,2.835e+00,2.960e+00,3.089e+00,
     &  3.221e+00,3.357e+00,3.498e+00,3.642e+00,3.790e+00,3.944e+00,
     &  4.104e+00,4.273e+00,4.445e+00,4.615e+00,4.791e+00,4.983e+00,
     &  5.195e+00,5.421e+00,5.656e+00,5.893e+00,6.127e+00,6.356e+00,
     &  6.582e+00,6.808e+00,7.036e+00,7.270e+00,7.517e+00,7.789e+00/
      data sun2a(121:180)/
     &  8.091e+00,8.407e+00,8.712e+00,8.990e+00,9.249e+00,9.500e+00,
     &  9.755e+00,1.001e+01,1.025e+01,1.048e+01,1.070e+01,1.095e+01,
     &  1.123e+01,1.155e+01,1.190e+01,1.225e+01,1.260e+01,1.293e+01,
     &  1.325e+01,1.353e+01,1.378e+01,1.404e+01,1.432e+01,1.466e+01,
     &  1.507e+01,1.553e+01,1.601e+01,1.643e+01,1.677e+01,1.708e+01,
     &  1.747e+01,1.796e+01,1.843e+01,1.873e+01,1.891e+01,1.914e+01,
     &  1.949e+01,1.984e+01,2.016e+01,2.051e+01,2.102e+01,2.177e+01,
     &  2.257e+01,2.319e+01,2.358e+01,2.390e+01,2.433e+01,2.483e+01,
     &  2.524e+01,2.565e+01,2.631e+01,2.721e+01,2.798e+01,2.842e+01,
     &  2.882e+01,2.957e+01,3.053e+01,3.125e+01,3.167e+01,3.222e+01/
      data sun2a(181:240)/
     &  3.309e+01,3.397e+01,3.460e+01,3.500e+01,3.540e+01,3.603e+01,
     &  3.699e+01,3.789e+01,3.840e+01,3.889e+01,3.986e+01,4.093e+01,
     &  4.157e+01,4.213e+01,4.308e+01,4.435e+01,4.552e+01,4.598e+01,
     &  4.628e+01,4.833e+01,5.199e+01,5.437e+01,5.408e+01,5.217e+01,
     &  5.071e+01,5.215e+01,5.571e+01,5.655e+01,5.441e+01,5.327e+01,
     &  5.608e+01,6.197e+01,6.441e+01,6.065e+01,5.515e+01,5.307e+01,
     &  5.748e+01,6.464e+01,6.835e+01,6.906e+01,6.987e+01,7.094e+01,
     &  7.166e+01,7.277e+01,7.433e+01,7.526e+01,7.488e+01,7.361e+01,
     &  7.321e+01,7.489e+01,7.804e+01,8.020e+01,8.088e+01,8.267e+01,
     &  8.498e+01,8.624e+01,8.836e+01,9.200e+01,9.538e+01,9.812e+01/
      data sun2a(241:300)/
     &  1.003e+02,1.006e+02,1.000e+02,1.018e+02,1.051e+02,1.075e+02,
     &  1.100e+02,1.124e+02,1.139e+02,1.138e+02,1.192e+02,1.220e+02,
     &  1.246e+02,1.271e+02,1.252e+02,1.244e+02,1.250e+02,1.279e+02,
     &  1.307e+02,1.320e+02,1.337e+02,1.367e+02,1.362e+02,1.350e+02,
     &  1.374e+02,1.384e+02,1.373e+02,1.364e+02,1.426e+02,1.445e+02,
     &  1.484e+02,1.519e+02,1.516e+02,1.554e+02,1.576e+02,1.597e+02,
     &  1.623e+02,1.684e+02,1.714e+02,1.698e+02,1.703e+02,1.723e+02,
     &  1.767e+02,1.819e+02,1.861e+02,1.879e+02,1.860e+02,1.898e+02,
     &  1.894e+02,1.929e+02,2.020e+02,2.096e+02,2.058e+02,2.129e+02,
     &  2.156e+02,2.165e+02,2.192e+02,2.203e+02,2.211e+02,2.271e+02/
      data sun2a(301:360)/
     &  2.300e+02,2.332e+02,2.339e+02,2.345e+02,2.344e+02,2.358e+02,
     &  2.398e+02,2.431e+02,2.412e+02,2.423e+02,2.437e+02,2.428e+02,
     &  2.462e+02,2.461e+02,2.468e+02,2.518e+02,2.554e+02,2.587e+02,
     &  2.603e+02,2.634e+02,2.687e+02,2.718e+02,2.730e+02,2.739e+02,
     &  2.747e+02,2.744e+02,2.797e+02,2.878e+02,2.877e+02,2.880e+02,
     &  2.900e+02,2.919e+02,2.953e+02,2.968e+02,3.005e+02,3.022e+02,
     &  2.991e+02,3.014e+02,3.057e+02,3.093e+02,3.106e+02,3.132e+02,
     &  3.146e+02,3.096e+02,3.188e+02,3.205e+02,3.216e+02,3.286e+02,
     &  3.317e+02,3.372e+02,3.456e+02,3.455e+02,3.430e+02,3.444e+02,
     &  3.462e+02,3.492e+02,3.518e+02,3.547e+02,3.570e+02,3.583e+02/
      data sun2a(361:420)/
     &  3.623e+02,3.641e+02,3.650e+02,3.678e+02,3.690e+02,3.691e+02,
     &  3.722e+02,3.778e+02,3.813e+02,3.842e+02,3.887e+02,3.936e+02,
     &  3.970e+02,3.987e+02,4.006e+02,4.041e+02,4.082e+02,4.125e+02,
     &  4.156e+02,4.162e+02,4.165e+02,4.195e+02,4.259e+02,4.333e+02,
     &  4.377e+02,4.381e+02,4.398e+02,4.415e+02,4.387e+02,4.343e+02,
     &  4.375e+02,4.490e+02,4.489e+02,4.395e+02,4.371e+02,4.393e+02,
     &  4.443e+02,4.550e+02,4.670e+02,4.730e+02,4.696e+02,4.675e+02,
     &  4.738e+02,4.775e+02,4.775e+02,4.810e+02,4.839e+02,4.822e+02,
     &  4.791e+02,4.821e+02,4.934e+02,4.984e+02,4.920e+02,4.895e+02,
     &  4.933e+02,4.955e+02,4.965e+02,4.996e+02,5.046e+02,5.097e+02/
      data sun2a(421:480)/
     &  5.120e+02,5.120e+02,5.123e+02,5.150e+02,5.207e+02,5.273e+02,
     &  5.319e+02,5.322e+02,5.305e+02,5.323e+02,5.393e+02,5.486e+02,
     &  5.530e+02,5.490e+02,5.460e+02,5.510e+02,5.564e+02,5.572e+02,
     &  5.578e+02,5.610e+02,5.640e+02,5.656e+02,5.664e+02,5.679e+02,
     &  5.715e+02,5.767e+02,5.815e+02,5.865e+02,5.936e+02,6.007e+02,
     &  6.028e+02,6.014e+02,6.030e+02,6.069e+02,6.060e+02,6.010e+02,
     &  6.008e+02,6.072e+02,6.129e+02,6.141e+02,6.144e+02,6.166e+02,
     &  6.205e+02,6.252e+02,6.298e+02,6.338e+02,6.373e+02,6.405e+02,
     &  6.425e+02,6.426e+02,6.419e+02,6.431e+02,6.467e+02,6.506e+02,
     &  6.543e+02,6.610e+02,6.721e+02,6.823e+02,6.849e+02,6.822e+02/
      data sun2a(481:540)/
     &  6.825e+02,6.878e+02,6.914e+02,6.896e+02,6.881e+02,6.937e+02,
     &  7.033e+02,7.081e+02,7.062e+02,7.046e+02,7.090e+02,7.173e+02,
     &  7.254e+02,7.311e+02,7.342e+02,7.354e+02,7.366e+02,7.393e+02,
     &  7.429e+02,7.450e+02,7.443e+02,7.424e+02,7.495e+02,7.557e+02,
     &  7.588e+02,7.663e+02,7.615e+02,7.621e+02,7.697e+02,7.642e+02,
     &  7.638e+02,7.689e+02,7.627e+02,7.539e+02,7.624e+02,7.658e+02,
     &  7.722e+02,7.607e+02,7.621e+02,7.668e+02,7.670e+02,7.693e+02,
     &  7.735e+02,7.668e+02,7.636e+02,7.738e+02,7.772e+02,7.796e+02,
     &  7.925e+02,7.975e+02,7.878e+02,7.938e+02,8.060e+02,8.048e+02,
     &  8.066e+02,8.217e+02,8.303e+02,8.275e+02,8.311e+02,8.302e+02/
      data sun2a(541:600)/
     &  8.262e+02,8.233e+02,8.222e+02,8.339e+02,8.546e+02,8.598e+02,
     &  8.626e+02,8.712e+02,8.752e+02,8.677e+02,8.639e+02,8.833e+02,
     &  8.934e+02,8.977e+02,9.052e+02,9.054e+02,9.111e+02,9.302e+02,
     &  9.392e+02,9.347e+02,9.352e+02,9.424e+02,9.481e+02,9.470e+02,
     &  9.519e+02,9.601e+02,9.519e+02,9.542e+02,9.591e+02,9.634e+02,
     &  9.802e+02,9.837e+02,9.788e+02,9.794e+02,9.852e+02,9.771e+02,
     &  9.199e+02,8.997e+02,9.629e+02,9.972e+02,9.999e+02,9.957e+02,
     &  9.999e+02,1.015e+03,9.516e+02,8.935e+02,9.551e+02,1.003e+03,
     &  9.901e+02,9.788e+02,1.011e+03,1.035e+03,1.032e+03,1.030e+03,
     &  1.040e+03,1.046e+03,1.044e+03,1.050e+03,1.056e+03,1.050e+03/
      data sun2a(601:660)/
     &  1.038e+03,1.052e+03,1.072e+03,1.076e+03,1.077e+03,1.079e+03,
     &  1.078e+03,1.076e+03,1.080e+03,1.081e+03,1.070e+03,1.078e+03,
     &  1.104e+03,1.111e+03,1.112e+03,1.118e+03,1.120e+03,1.109e+03,
     &  1.101e+03,1.113e+03,1.123e+03,1.120e+03,1.124e+03,1.136e+03,
     &  1.144e+03,1.141e+03,1.141e+03,1.152e+03,1.149e+03,1.138e+03,
     &  1.141e+03,1.151e+03,1.160e+03,1.171e+03,1.178e+03,1.180e+03,
     &  1.182e+03,1.182e+03,1.180e+03,1.182e+03,1.188e+03,1.190e+03,
     &  1.191e+03,1.197e+03,1.196e+03,1.192e+03,1.201e+03,1.210e+03,
     &  1.209e+03,1.208e+03,1.205e+03,1.193e+03,1.193e+03,1.220e+03,
     &  1.243e+03,1.245e+03,1.242e+03,1.240e+03,1.241e+03,1.244e+03/
      data sun2a(661:720)/
     &  1.249e+03,1.253e+03,1.257e+03,1.260e+03,1.262e+03,1.264e+03,
     &  1.266e+03,1.270e+03,1.277e+03,1.284e+03,1.284e+03,1.283e+03,
     &  1.287e+03,1.287e+03,1.272e+03,1.262e+03,1.271e+03,1.289e+03,
     &  1.305e+03,1.312e+03,1.312e+03,1.314e+03,1.320e+03,1.326e+03,
     &  1.328e+03,1.325e+03,1.323e+03,1.325e+03,1.335e+03,1.346e+03,
     &  1.354e+03,1.354e+03,1.347e+03,1.338e+03,1.331e+03,1.330e+03,
     &  1.338e+03,1.352e+03,1.363e+03,1.369e+03,1.372e+03,1.376e+03,
     &  1.382e+03,1.388e+03,1.389e+03,1.388e+03,1.392e+03,1.402e+03,
     &  1.413e+03,1.418e+03,1.411e+03,1.396e+03,1.386e+03,1.388e+03,
     &  1.405e+03,1.424e+03,1.428e+03,1.422e+03,1.424e+03,1.435e+03/
      data sun2a(721:780)/
     &  1.445e+03,1.451e+03,1.452e+03,1.452e+03,1.454e+03,1.460e+03,
     &  1.467e+03,1.471e+03,1.469e+03,1.463e+03,1.460e+03,1.469e+03,
     &  1.482e+03,1.491e+03,1.495e+03,1.498e+03,1.501e+03,1.505e+03,
     &  1.510e+03,1.512e+03,1.513e+03,1.516e+03,1.522e+03,1.524e+03,
     &  1.521e+03,1.520e+03,1.532e+03,1.546e+03,1.548e+03,1.542e+03,
     &  1.542e+03,1.554e+03,1.564e+03,1.564e+03,1.560e+03,1.561e+03,
     &  1.570e+03,1.582e+03,1.578e+03,1.530e+03,1.447e+03,1.397e+03,
     &  1.429e+03,1.506e+03,1.567e+03,1.594e+03,1.606e+03,1.614e+03,
     &  1.609e+03,1.589e+03,1.568e+03,1.567e+03,1.587e+03,1.610e+03,
     &  1.624e+03,1.630e+03,1.631e+03,1.628e+03,1.622e+03,1.617e+03/
      data sun2a(781:840)/
     &  1.619e+03,1.632e+03,1.648e+03,1.658e+03,1.660e+03,1.658e+03,
     &  1.658e+03,1.659e+03,1.660e+03,1.659e+03,1.654e+03,1.645e+03,
     &  1.642e+03,1.653e+03,1.674e+03,1.694e+03,1.701e+03,1.703e+03,
     &  1.698e+03,1.655e+03,1.644e+03,1.662e+03,1.676e+03,1.708e+03,
     &  1.703e+03,1.711e+03,1.732e+03,1.717e+03,1.720e+03,1.730e+03,
     &  1.683e+03,1.629e+03,1.684e+03,1.727e+03,1.708e+03,1.689e+03,
     &  1.698e+03,1.733e+03,1.738e+03,1.714e+03,1.735e+03,1.750e+03,
     &  1.750e+03,1.760e+03,1.764e+03,1.765e+03,1.769e+03,1.780e+03,
     &  1.793e+03,1.765e+03,1.729e+03,1.746e+03,1.753e+03,1.758e+03,
     &  1.775e+03,1.768e+03,1.768e+03,1.790e+03,1.807e+03,1.799e+03/
      data sun2a(841:900)/
     &  1.783e+03,1.779e+03,1.792e+03,1.810e+03,1.808e+03,1.794e+03,
     &  1.819e+03,1.774e+03,1.649e+03,1.674e+03,1.789e+03,1.847e+03,
     &  1.848e+03,1.813e+03,1.796e+03,1.840e+03,1.868e+03,1.865e+03,
     &  1.873e+03,1.872e+03,1.856e+03,1.845e+03,1.842e+03,1.824e+03,
     &  1.795e+03,1.820e+03,1.862e+03,1.858e+03,1.839e+03,1.841e+03,
     &  1.864e+03,1.877e+03,1.884e+03,1.895e+03,1.875e+03,1.821e+03,
     &  1.779e+03,1.810e+03,1.855e+03,1.832e+03,1.837e+03,1.882e+03,
     &  1.866e+03,1.820e+03,1.805e+03,1.831e+03,1.862e+03,1.867e+03,
     &  1.863e+03,1.852e+03,1.835e+03,1.835e+03,1.845e+03,1.832e+03,
     &  1.804e+03,1.793e+03,1.822e+03,1.846e+03,1.832e+03,1.848e+03/
      data sun2a(901:960)/
     &  1.894e+03,1.909e+03,1.901e+03,1.891e+03,1.870e+03,1.854e+03,
     &  1.866e+03,1.874e+03,1.869e+03,1.882e+03,1.897e+03,1.884e+03,
     &  1.856e+03,1.841e+03,1.855e+03,1.885e+03,1.904e+03,1.900e+03,
     &  1.887e+03,1.888e+03,1.879e+03,1.845e+03,1.844e+03,1.877e+03,
     &  1.847e+03,1.785e+03,1.793e+03,1.849e+03,1.894e+03,1.909e+03,
     &  1.893e+03,1.867e+03,1.886e+03,1.960e+03,1.972e+03,1.896e+03,
     &  1.884e+03,1.918e+03,1.854e+03,1.793e+03,1.876e+03,1.974e+03,
     &  1.976e+03,1.944e+03,1.926e+03,1.914e+03,1.903e+03,1.883e+03,
     &  1.813e+03,1.711e+03,1.718e+03,1.860e+03,1.965e+03,1.970e+03,
     &  1.941e+03,1.903e+03,1.852e+03,1.836e+03,1.879e+03,1.902e+03/
      data sun2a(961:1020)/
     &  1.863e+03,1.839e+03,1.841e+03,1.780e+03,1.685e+03,1.677e+03,
     &  1.719e+03,1.697e+03,1.684e+03,1.785e+03,1.898e+03,1.910e+03,
     &  1.877e+03,1.867e+03,1.863e+03,1.860e+03,1.900e+03,1.971e+03,
     &  2.000e+03,1.971e+03,1.937e+03,1.923e+03,1.923e+03,1.924e+03,
     &  1.917e+03,1.912e+03,1.926e+03,1.960e+03,1.995e+03,1.996e+03,
     &  1.939e+03,1.884e+03,1.895e+03,1.933e+03,1.935e+03,1.899e+03,
     &  1.853e+03,1.820e+03,1.822e+03,1.865e+03,1.936e+03,1.966e+03,
     &  1.920e+03,1.881e+03,1.932e+03,2.016e+03,2.050e+03,2.021e+03,
     &  1.961e+03,1.938e+03,1.997e+03,2.051e+03,2.003e+03,1.912e+03,
     &  1.880e+03,1.895e+03,1.898e+03,1.899e+03,1.938e+03,1.994e+03/
      data sun2a(1021:1080)/
     &  2.010e+03,1.982e+03,1.949e+03,1.927e+03,1.912e+03,1.878e+03,
     &  1.792e+03,1.680e+03,1.645e+03,1.727e+03,1.845e+03,1.926e+03,
     &  1.973e+03,2.005e+03,2.022e+03,2.022e+03,2.026e+03,2.054e+03,
     &  2.087e+03,2.083e+03,2.053e+03,2.047e+03,2.070e+03,2.072e+03,
     &  2.038e+03,2.020e+03,2.050e+03,2.074e+03,2.038e+03,1.979e+03,
     &  1.964e+03,1.997e+03,2.038e+03,2.058e+03,2.048e+03,2.018e+03,
     &  1.999e+03,2.011e+03,2.040e+03,2.056e+03,2.040e+03,1.982e+03,
     &  1.911e+03,1.892e+03,1.938e+03,1.992e+03,2.006e+03,2.001e+03,
     &  2.011e+03,2.023e+03,1.998e+03,1.948e+03,1.936e+03,1.987e+03,
     &  2.038e+03,2.033e+03,1.996e+03,1.984e+03,2.012e+03,2.056e+03/
      data sun2a(1081:1140)/
     &  2.092e+03,2.107e+03,2.095e+03,2.070e+03,2.053e+03,2.047e+03,
     &  2.044e+03,2.036e+03,2.017e+03,1.988e+03,1.973e+03,1.999e+03,
     &  2.057e+03,2.104e+03,2.109e+03,2.089e+03,2.069e+03,2.052e+03,
     &  2.031e+03,2.006e+03,1.987e+03,1.982e+03,1.979e+03,1.964e+03,
     &  1.944e+03,1.952e+03,2.007e+03,2.083e+03,2.139e+03,2.158e+03,
     &  2.143e+03,2.103e+03,2.051e+03,2.002e+03,1.975e+03,1.988e+03,
     &  2.038e+03,2.075e+03,2.051e+03,1.972e+03,1.885e+03,1.829e+03,
     &  1.821e+03,1.866e+03,1.935e+03,1.974e+03,1.959e+03,1.925e+03,
     &  1.920e+03,1.950e+03,1.985e+03,1.996e+03,1.966e+03,1.885e+03,
     &  1.782e+03,1.727e+03,1.759e+03,1.817e+03,1.800e+03,1.693e+03/
      data sun2a(1141:1200)/
     &  1.593e+03,1.599e+03,1.700e+03,1.824e+03,1.910e+03,1.938e+03,
     &  1.903e+03,1.822e+03,1.738e+03,1.683e+03,1.667e+03,1.683e+03,
     &  1.715e+03,1.734e+03,1.712e+03,1.668e+03,1.655e+03,1.698e+03,
     &  1.727e+03,1.637e+03,1.416e+03,1.204e+03,1.156e+03,1.278e+03,
     &  1.450e+03,1.561e+03,1.595e+03,1.588e+03,1.571e+03,1.566e+03,
     &  1.590e+03,1.641e+03,1.688e+03,1.708e+03,1.704e+03,1.701e+03,
     &  1.719e+03,1.749e+03,1.772e+03,1.773e+03,1.745e+03,1.690e+03,
     &  1.625e+03,1.589e+03,1.619e+03,1.701e+03,1.783e+03,1.816e+03,
     &  1.801e+03,1.765e+03,1.734e+03,1.715e+03,1.705e+03,1.702e+03,
     &  1.697e+03,1.682e+03,1.661e+03,1.657e+03,1.693e+03,1.763e+03/
      data sun2a(1201:1260)/
     &  1.827e+03,1.842e+03,1.806e+03,1.756e+03,1.726e+03,1.724e+03,
     &  1.737e+03,1.749e+03,1.756e+03,1.760e+03,1.762e+03,1.770e+03,
     &  1.792e+03,1.827e+03,1.849e+03,1.820e+03,1.721e+03,1.596e+03,
     &  1.514e+03,1.523e+03,1.602e+03,1.706e+03,1.793e+03,1.838e+03,
     &  1.820e+03,1.738e+03,1.631e+03,1.553e+03,1.539e+03,1.574e+03,
     &  1.624e+03,1.661e+03,1.677e+03,1.673e+03,1.653e+03,1.626e+03,
     &  1.607e+03,1.604e+03,1.621e+03,1.655e+03,1.701e+03,1.752e+03,
     &  1.796e+03,1.823e+03,1.827e+03,1.809e+03,1.767e+03,1.714e+03,
     &  1.667e+03,1.644e+03,1.644e+03,1.653e+03,1.655e+03,1.639e+03,
     &  1.592e+03,1.506e+03,1.377e+03,1.210e+03,1.011e+03,8.076e+02/
      data sun2a(1261:1320)/
     &  6.668e+02,6.645e+02,8.352e+02,1.100e+03,1.331e+03,1.423e+03,
     &  1.364e+03,1.194e+03,9.618e+02,7.250e+02,5.513e+02,5.040e+02,
     &  5.963e+02,7.752e+02,9.756e+02,1.150e+03,1.287e+03,1.386e+03,
     &  1.448e+03,1.474e+03,1.469e+03,1.435e+03,1.377e+03,1.296e+03,
     &  1.196e+03,1.085e+03,9.854e+02,9.173e+02,8.946e+02,9.109e+02,
     &  9.515e+02,1.002e+03,1.046e+03,1.071e+03,1.061e+03,1.021e+03,
     &  9.772e+02,9.592e+02,9.821e+02,1.021e+03,1.033e+03,9.834e+02,
     &  8.798e+02,7.627e+02,6.753e+02,6.433e+02,6.627e+02,7.215e+02,
     &  8.083e+02,9.132e+02,1.027e+03,1.140e+03,1.236e+03,1.293e+03,
     &  1.287e+03,1.210e+03,1.102e+03,1.022e+03,1.023e+03,1.109e+03/
      data sun2a(1321:1380)/
     &  1.233e+03,1.337e+03,1.383e+03,1.373e+03,1.325e+03,1.258e+03,
     &  1.189e+03,1.134e+03,1.107e+03,1.114e+03,1.137e+03,1.148e+03,
     &  1.121e+03,1.054e+03,9.681e+02,8.892e+02,8.379e+02,8.176e+02,
     &  8.237e+02,8.510e+02,8.965e+02,9.598e+02,1.041e+03,1.138e+03,
     &  1.231e+03,1.294e+03,1.300e+03,1.241e+03,1.155e+03,1.092e+03,
     &  1.097e+03,1.170e+03,1.264e+03,1.322e+03,1.307e+03,1.234e+03,
     &  1.146e+03,1.091e+03,1.093e+03,1.135e+03,1.189e+03,1.229e+03,
     &  1.246e+03,1.249e+03,1.250e+03,1.261e+03,1.275e+03,1.280e+03,
     &  1.262e+03,1.214e+03,1.145e+03,1.070e+03,1.001e+03,9.525e+02,
     &  9.305e+02,9.417e+02,9.903e+02,1.064e+03,1.135e+03,1.172e+03/
      data sun2a(1381:1440)/
     &  1.149e+03,1.076e+03,9.843e+02,9.063e+02,8.682e+02,8.738e+02,
     &  9.153e+02,9.844e+02,1.067e+03,1.137e+03,1.163e+03,1.116e+03,
     &  9.905e+02,8.309e+02,6.923e+02,6.274e+02,6.541e+02,7.392e+02,
     &  8.389e+02,9.117e+02,9.419e+02,9.444e+02,9.396e+02,9.461e+02,
     &  9.702e+02,1.005e+03,1.042e+03,1.074e+03,1.097e+03,1.114e+03,
     &  1.129e+03,1.143e+03,1.153e+03,1.152e+03,1.132e+03,1.084e+03,
     &  1.017e+03,9.460e+02,8.904e+02,8.662e+02,8.765e+02,9.131e+02,
     &  9.661e+02,1.025e+03,1.080e+03,1.119e+03,1.103e+03,1.244e+03,
     &  1.210e+03,1.079e+03,8.522e+02,9.568e+02,8.423e+02,8.974e+02,
     &  1.082e+03,9.142e+02,9.931e+02,1.050e+03,8.450e+02,8.392e+02/

      data (sun2b(i),i=1,60)/
     &  8.765e+02,8.922e+02,9.131e+02,9.382e+02,9.661e+02,9.956e+02,
     &  1.025e+03,1.054e+03,1.080e+03,1.102e+03,1.119e+03,1.132e+03,
     &  1.103e+03,1.159e+03,1.244e+03,1.238e+03,1.210e+03,1.196e+03,
     &  1.079e+03,8.956e+02,8.522e+02,9.356e+02,9.568e+02,8.971e+02,
     &  8.423e+02,8.212e+02,8.974e+02,1.043e+03,1.082e+03,9.888e+02,
     &  9.142e+02,9.294e+02,9.931e+02,1.042e+03,1.050e+03,9.843e+02,
     &  8.450e+02,7.708e+02,8.392e+02,9.397e+02,1.026e+03,1.121e+03,
     &  1.163e+03,1.143e+03,1.078e+03,1.027e+03,1.078e+03,1.094e+03,
     &  9.698e+02,8.537e+02,8.499e+02,9.091e+02,9.957e+02,1.095e+03,
     &  1.147e+03,1.086e+03,1.010e+03,1.065e+03,1.129e+03,1.081e+03/
      data (sun2b(i),i=61,120)/
     &  9.879e+02,8.982e+02,8.352e+02,7.716e+02,6.871e+02,6.145e+02,
     &  6.061e+02,7.371e+02,9.081e+02,9.976e+02,1.081e+03,1.126e+03,
     &  1.057e+03,1.028e+03,1.142e+03,1.253e+03,1.225e+03,1.103e+03,
     &  1.039e+03,1.043e+03,1.003e+03,9.655e+02,1.035e+03,1.151e+03,
     &  1.201e+03,1.152e+03,1.069e+03,9.958e+02,8.895e+02,8.185e+02,
     &  9.070e+02,1.042e+03,1.056e+03,1.001e+03,9.720e+02,9.857e+02,
     &  1.027e+03,1.055e+03,1.078e+03,1.127e+03,1.205e+03,1.246e+03,
     &  1.201e+03,1.145e+03,1.098e+03,1.030e+03,9.268e+02,8.367e+02,
     &  8.641e+02,9.935e+02,1.075e+03,1.033e+03,1.009e+03,1.066e+03,
     &  1.067e+03,1.005e+03,9.715e+02,9.232e+02,8.157e+02,7.997e+02/
      data (sun2b(i),i=121,180)/
     &  9.462e+02,1.100e+03,1.126e+03,1.032e+03,8.951e+02,7.843e+02,
     &  7.348e+02,7.265e+02,7.269e+02,7.655e+02,8.639e+02,9.922e+02,
     &  1.071e+03,1.028e+03,8.588e+02,6.472e+02,5.632e+02,6.800e+02,
     &  9.064e+02,1.094e+03,1.155e+03,1.124e+03,1.098e+03,1.110e+03,
     &  1.076e+03,9.442e+02,8.492e+02,9.285e+02,1.062e+03,1.119e+03,
     &  1.119e+03,1.075e+03,1.006e+03,9.800e+02,9.991e+02,1.002e+03,
     &  9.398e+02,8.381e+02,8.161e+02,9.087e+02,1.015e+03,1.058e+03,
     &  1.044e+03,9.875e+02,9.463e+02,9.814e+02,1.056e+03,1.094e+03,
     &  1.028e+03,9.164e+02,9.090e+02,9.918e+02,1.050e+03,1.076e+03,
     &  1.094e+03,1.076e+03,1.015e+03,9.496e+02,9.473e+02,1.001e+03/
      data (sun2b(i),i=181,240)/
     &  1.052e+03,1.073e+03,1.068e+03,1.013e+03,9.078e+02,8.663e+02,
     &  9.509e+02,1.038e+03,1.080e+03,1.184e+03,1.291e+03,1.269e+03,
     &  1.199e+03,1.189e+03,1.188e+03,1.187e+03,1.198e+03,1.171e+03,
     &  1.133e+03,1.132e+03,1.096e+03,9.711e+02,8.471e+02,8.366e+02,
     &  9.228e+02,9.910e+02,9.875e+02,9.692e+02,9.815e+02,9.814e+02,
     &  9.720e+02,9.853e+02,1.003e+03,1.037e+03,1.071e+03,1.066e+03,
     &  1.027e+03,9.848e+02,1.003e+03,1.070e+03,1.118e+03,1.116e+03,
     &  1.049e+03,9.653e+02,9.723e+02,1.046e+03,1.097e+03,1.128e+03,
     &  1.134e+03,1.100e+03,1.079e+03,1.083e+03,1.027e+03,9.275e+02,
     &  8.791e+02,8.588e+02,8.310e+02,8.078e+02,7.896e+02,8.138e+02/
      data (sun2b(i),i=241,300)/
     &  8.935e+02,9.376e+02,9.016e+02,8.645e+02,8.733e+02,8.910e+02,
     &  8.625e+02,8.103e+02,7.874e+02,7.529e+02,7.153e+02,7.081e+02,
     &  7.289e+02,7.868e+02,8.077e+02,7.363e+02,6.451e+02,6.169e+02,
     &  6.492e+02,6.918e+02,7.492e+02,8.202e+02,8.207e+02,7.913e+02,
     &  8.543e+02,9.406e+02,9.564e+02,9.094e+02,8.242e+02,7.672e+02,
     &  7.221e+02,6.534e+02,6.247e+02,6.337e+02,6.551e+02,7.079e+02,
     &  7.849e+02,8.808e+02,9.612e+02,9.856e+02,9.862e+02,9.665e+02,
     &  9.215e+02,8.889e+02,8.558e+02,8.517e+02,8.868e+02,8.510e+02,
     &  7.670e+02,7.390e+02,7.245e+02,6.576e+02,5.878e+02,6.169e+02,
     &  7.606e+02,9.032e+02,9.173e+02,8.385e+02,7.848e+02,7.594e+02/
      data (sun2b(i),i=301,360)/
     &  7.196e+02,6.715e+02,6.246e+02,5.886e+02,5.747e+02,5.967e+02,
     &  6.980e+02,8.664e+02,9.748e+02,9.604e+02,9.301e+02,9.627e+02,
     &  1.007e+03,1.002e+03,9.263e+02,8.166e+02,7.633e+02,7.729e+02,
     &  7.627e+02,7.294e+02,7.250e+02,7.272e+02,6.727e+02,5.814e+02,
     &  5.210e+02,4.888e+02,4.786e+02,5.421e+02,6.637e+02,7.495e+02,
     &  7.859e+02,8.110e+02,8.182e+02,8.138e+02,8.245e+02,8.366e+02,
     &  7.997e+02,7.280e+02,6.604e+02,5.593e+02,4.733e+02,5.502e+02,
     &  7.520e+02,8.858e+02,9.068e+02,9.122e+02,9.293e+02,8.997e+02,
     &  8.302e+02,7.746e+02,7.364e+02,7.241e+02,7.401e+02,7.541e+02,
     &  7.650e+02,7.808e+02,7.889e+02,7.849e+02,7.588e+02,7.259e+02/
      data (sun2b(i),i=361,420)/
     &  7.518e+02,8.042e+02,7.777e+02,7.034e+02,6.653e+02,6.640e+02,
     &  6.794e+02,7.061e+02,7.576e+02,8.361e+02,8.800e+02,8.812e+02,
     &  9.079e+02,9.293e+02,8.943e+02,8.740e+02,9.186e+02,9.535e+02,
     &  9.223e+02,8.666e+02,8.365e+02,8.253e+02,7.525e+02,5.860e+02,
     &  4.275e+02,3.740e+02,4.372e+02,5.343e+02,5.567e+02,5.631e+02,
     &  6.293e+02,6.313e+02,5.188e+02,4.383e+02,4.603e+02,5.305e+02,
     &  6.085e+02,6.580e+02,6.621e+02,6.862e+02,7.752e+02,8.431e+02,
     &  7.975e+02,6.853e+02,6.113e+02,6.287e+02,7.114e+02,7.549e+02,
     &  7.288e+02,7.228e+02,7.264e+02,6.797e+02,6.658e+02,7.105e+02,
     &  7.231e+02,7.241e+02,7.602e+02,7.840e+02,7.428e+02,6.343e+02/
      data (sun2b(i),i=421,480)/
     &  5.465e+02,5.635e+02,6.110e+02,6.232e+02,6.654e+02,7.435e+02,
     &  7.645e+02,6.711e+02,5.132e+02,4.019e+02,4.058e+02,5.157e+02,
     &  6.399e+02,6.778e+02,6.795e+02,7.593e+02,8.481e+02,8.199e+02,
     &  7.518e+02,7.105e+02,6.153e+02,5.251e+02,5.833e+02,7.152e+02,
     &  7.675e+02,7.391e+02,6.640e+02,5.806e+02,5.728e+02,6.341e+02,
     &  6.488e+02,5.613e+02,4.977e+02,5.917e+02,7.378e+02,7.942e+02,
     &  8.025e+02,7.993e+02,7.358e+02,6.584e+02,6.595e+02,7.182e+02,
     &  7.617e+02,6.972e+02,5.451e+02,4.745e+02,5.270e+02,5.977e+02,
     &  5.847e+02,4.473e+02,2.914e+02,2.613e+02,3.303e+02,4.020e+02,
     &  4.663e+02,5.313e+02,5.723e+02,5.849e+02,5.852e+02,5.695e+02/
      data (sun2b(i),i=481,540)/
     &  5.583e+02,5.594e+02,5.120e+02,4.264e+02,3.781e+02,3.983e+02,
     &  4.735e+02,5.422e+02,5.318e+02,4.375e+02,3.419e+02,3.058e+02,
     &  2.999e+02,3.281e+02,4.400e+02,5.865e+02,6.603e+02,6.252e+02,
     &  5.103e+02,4.189e+02,4.474e+02,5.349e+02,6.059e+02,6.671e+02,
     &  6.873e+02,6.368e+02,5.496e+02,4.729e+02,4.195e+02,3.701e+02,
     &  3.280e+02,3.205e+02,3.540e+02,3.992e+02,4.510e+02,5.283e+02,
     &  6.083e+02,6.961e+02,7.743e+02,7.608e+02,6.906e+02,6.482e+02,
     &  5.806e+02,4.780e+02,4.539e+02,4.887e+02,4.640e+02,4.216e+02,
     &  4.443e+02,4.466e+02,3.760e+02,3.421e+02,3.975e+02,5.110e+02,
     &  6.464e+02,7.251e+02,7.031e+02,6.391e+02,6.191e+02,6.547e+02/
      data (sun2b(i),i=541,600)/
     &  6.660e+02,6.114e+02,5.802e+02,6.073e+02,5.910e+02,5.423e+02,
     &  5.838e+02,6.730e+02,6.732e+02,5.824e+02,4.657e+02,3.773e+02,
     &  3.770e+02,4.873e+02,6.079e+02,6.175e+02,5.835e+02,6.017e+02,
     &  6.159e+02,5.755e+02,5.416e+02,5.421e+02,5.223e+02,4.725e+02,
     &  4.233e+02,4.381e+02,5.567e+02,6.643e+02,6.699e+02,6.575e+02,
     &  6.847e+02,7.057e+02,6.831e+02,6.008e+02,5.099e+02,4.976e+02,
     &  5.111e+02,4.961e+02,5.003e+02,5.187e+02,5.299e+02,5.630e+02,
     &  6.092e+02,6.265e+02,6.221e+02,6.157e+02,6.004e+02,5.913e+02,
     &  5.981e+02,5.931e+02,5.909e+02,6.316e+02,6.965e+02,7.185e+02,
     &  6.761e+02,6.316e+02,6.196e+02,6.205e+02,6.241e+02,6.366e+02/
      data (sun2b(i),i=601,660)/
     &  6.580e+02,6.888e+02,7.248e+02,7.426e+02,7.223e+02,6.759e+02,
     &  6.660e+02,7.047e+02,7.037e+02,6.450e+02,5.983e+02,5.878e+02,
     &  5.909e+02,5.759e+02,5.280e+02,4.779e+02,4.575e+02,4.568e+02,
     &  4.549e+02,4.486e+02,4.455e+02,4.454e+02,4.444e+02,4.460e+02,
     &  4.559e+02,4.680e+02,4.543e+02,3.933e+02,3.012e+02,2.114e+02,
     &  1.671e+02,1.940e+02,2.540e+02,3.054e+02,3.530e+02,3.851e+02,
     &  3.870e+02,3.916e+02,4.062e+02,4.153e+02,4.353e+02,4.698e+02,
     &  4.921e+02,4.727e+02,4.099e+02,3.533e+02,3.407e+02,3.553e+02,
     &  3.798e+02,4.018e+02,4.097e+02,4.069e+02,3.932e+02,3.789e+02,
     &  3.752e+02,3.735e+02,3.602e+02,3.228e+02,2.735e+02,2.378e+02/
      data (sun2b(i),i=661,720)/
     &  2.123e+02,1.848e+02,1.562e+02,1.278e+02,9.627e+01,6.881e+01,
     &  6.205e+01,7.714e+01,1.005e+02,1.276e+02,1.599e+02,1.941e+02,
     &  2.252e+02,2.546e+02,2.858e+02,3.001e+02,2.944e+02,3.089e+02,
     &  3.408e+02,3.463e+02,3.363e+02,3.475e+02,3.738e+02,3.888e+02,
     &  3.727e+02,3.253e+02,2.944e+02,3.176e+02,3.603e+02,3.781e+02,
     &  3.742e+02,3.740e+02,3.833e+02,3.879e+02,3.775e+02,3.570e+02,
     &  3.407e+02,3.287e+02,3.140e+02,3.169e+02,3.445e+02,3.555e+02,
     &  3.357e+02,3.187e+02,3.186e+02,3.224e+02,3.186e+02,3.049e+02,
     &  2.848e+02,2.681e+02,2.658e+02,2.735e+02,2.742e+02,2.524e+02,
     &  2.150e+02,1.886e+02,1.813e+02,1.813e+02,1.808e+02,1.752e+02/
      data (sun2b(i),i=721,780)/
     &  1.621e+02,1.451e+02,1.288e+02,1.138e+02,9.808e+01,8.307e+01,
     &  7.622e+01,7.836e+01,7.843e+01,7.424e+01,7.584e+01,8.032e+01,
     &  7.786e+01,7.030e+01,6.465e+01,6.705e+01,7.781e+01,8.317e+01,
     &  7.529e+01,7.120e+01,8.055e+01,9.201e+01,1.002e+02,1.086e+02,
     &  1.194e+02,1.308e+02,1.423e+02,1.589e+02,1.771e+02,1.864e+02,
     &  1.866e+02,1.815e+02,1.753e+02,1.755e+02,1.790e+02,1.770e+02,
     &  1.726e+02,1.727e+02,1.790e+02,1.938e+02,2.151e+02,2.336e+02,
     &  2.521e+02,2.777e+02,2.989e+02,2.984e+02,2.808e+02,2.742e+02,
     &  2.865e+02,2.855e+02,2.597e+02,2.414e+02,2.470e+02,2.599e+02,
     &  2.743e+02,2.985e+02,3.169e+02,3.032e+02,2.637e+02,2.293e+02/
      data (sun2b(i),i=781,840)/
     &  2.279e+02,2.561e+02,2.816e+02,3.002e+02,3.106e+02,2.795e+02,
     &  2.119e+02,1.522e+02,1.299e+02,1.475e+02,1.816e+02,2.154e+02,
     &  2.395e+02,2.331e+02,1.916e+02,1.394e+02,1.105e+02,1.189e+02,
     &  1.348e+02,1.291e+02,1.244e+02,1.435e+02,1.583e+02,1.418e+02,
     &  1.163e+02,1.116e+02,1.289e+02,1.492e+02,1.534e+02,1.456e+02,
     &  1.485e+02,1.593e+02,1.558e+02,1.542e+02,1.773e+02,2.034e+02,
     &  2.074e+02,2.053e+02,2.229e+02,2.532e+02,2.713e+02,2.793e+02,
     &  3.022e+02,3.215e+02,2.888e+02,2.301e+02,2.064e+02,2.132e+02,
     &  2.165e+02,2.075e+02,1.962e+02,1.952e+02,2.020e+02,1.943e+02,
     &  1.649e+02,1.366e+02,1.239e+02,1.281e+02,1.619e+02,2.170e+02/
      data (sun2b(i),i=841,900)/
     &  2.537e+02,2.493e+02,2.229e+02,2.131e+02,2.436e+02,2.931e+02,
     &  3.094e+02,2.864e+02,2.696e+02,2.722e+02,2.717e+02,2.658e+02,
     &  2.656e+02,2.648e+02,2.660e+02,2.895e+02,3.257e+02,3.373e+02,
     &  3.212e+02,3.003e+02,2.826e+02,2.871e+02,3.221e+02,3.358e+02,
     &  2.972e+02,2.541e+02,2.435e+02,2.395e+02,2.193e+02,2.119e+02,
     &  2.393e+02,2.714e+02,2.794e+02,2.723e+02,2.648e+02,2.505e+02,
     &  2.299e+02,2.221e+02,2.353e+02,2.568e+02,2.753e+02,2.869e+02,
     &  2.849e+02,2.695e+02,2.551e+02,2.535e+02,2.632e+02,2.748e+02,
     &  2.792e+02,2.702e+02,2.494e+02,2.290e+02,2.216e+02,2.314e+02,
     &  2.527e+02,2.806e+02,3.101e+02,3.283e+02,3.250e+02,2.903e+02/
      data (sun2b(i),i=901,960)/
     &  2.390e+02,2.234e+02,2.572e+02,2.826e+02,2.643e+02,2.433e+02,
     &  2.532e+02,2.729e+02,2.713e+02,2.561e+02,2.602e+02,2.714e+02,
     &  2.571e+02,2.366e+02,2.387e+02,2.489e+02,2.559e+02,2.720e+02,
     &  2.918e+02,2.974e+02,2.881e+02,2.833e+02,2.929e+02,3.017e+02,
     &  3.091e+02,3.220e+02,3.204e+02,2.954e+02,2.696e+02,2.544e+02,
     &  2.409e+02,2.282e+02,2.212e+02,2.137e+02,2.012e+02,1.972e+02,
     &  2.123e+02,2.334e+02,2.476e+02,2.617e+02,2.862e+02,3.225e+02,
     &  3.495e+02,3.383e+02,2.971e+02,2.615e+02,2.523e+02,2.646e+02,
     &  2.869e+02,2.989e+02,2.805e+02,2.444e+02,2.135e+02,1.930e+02,
     &  1.821e+02,1.685e+02,1.431e+02,1.141e+02,8.961e+01,7.359e+01/
      data (sun2b(i),i=961,1020)/
     &  7.399e+01,8.791e+01,9.626e+01,9.481e+01,9.660e+01,1.023e+02,
     &  1.022e+02,1.031e+02,1.178e+02,1.374e+02,1.461e+02,1.443e+02,
     &  1.379e+02,1.281e+02,1.228e+02,1.282e+02,1.307e+02,1.173e+02,
     &  9.891e+01,9.340e+01,1.056e+02,1.227e+02,1.264e+02,1.131e+02,
     &  9.232e+01,7.634e+01,6.903e+01,6.632e+01,7.128e+01,8.743e+01,
     &  1.059e+02,1.140e+02,1.079e+02,9.187e+01,7.521e+01,6.912e+01,
     &  7.593e+01,9.093e+01,1.097e+02,1.257e+02,1.358e+02,1.411e+02,
     &  1.381e+02,1.213e+02,9.181e+01,6.350e+01,5.211e+01,5.956e+01,
     &  8.101e+01,1.067e+02,1.190e+02,1.164e+02,1.108e+02,1.009e+02,
     &  8.906e+01,9.043e+01,1.044e+02,1.149e+02,1.248e+02,1.489e+02/
      data (sun2b(i),i=1021,1080)/
     &  1.717e+02,1.672e+02,1.423e+02,1.184e+02,9.865e+01,7.891e+01,
     &  6.813e+01,7.729e+01,1.009e+02,1.201e+02,1.255e+02,1.318e+02,
     &  1.557e+02,1.808e+02,1.818e+02,1.668e+02,1.501e+02,1.332e+02,
     &  1.161e+02,9.773e+01,8.163e+01,7.669e+01,8.761e+01,1.102e+02,
     &  1.349e+02,1.491e+02,1.476e+02,1.399e+02,1.352e+02,1.351e+02,
     &  1.380e+02,1.367e+02,1.288e+02,1.222e+02,1.205e+02,1.220e+02,
     &  1.231e+02,1.163e+02,1.014e+02,8.630e+01,7.472e+01,6.880e+01,
     &  7.133e+01,8.063e+01,9.049e+01,9.674e+01,1.007e+02,1.008e+02,
     &  9.368e+01,8.474e+01,8.153e+01,8.289e+01,8.456e+01,8.758e+01,
     &  9.178e+01,9.127e+01,8.701e+01,8.739e+01,9.015e+01,8.492e+01/
      data (sun2b(i),i=1081,1140)/
     &  7.127e+01,5.787e+01,5.186e+01,5.388e+01,5.791e+01,5.851e+01,
     &  5.702e+01,5.743e+01,6.067e+01,6.467e+01,6.736e+01,6.751e+01,
     &  6.423e+01,5.903e+01,5.570e+01,5.664e+01,5.940e+01,5.907e+01,
     &  5.652e+01,5.583e+01,5.586e+01,5.404e+01,5.198e+01,5.234e+01,
     &  5.467e+01,5.645e+01,5.675e+01,5.677e+01,5.800e+01,6.003e+01,
     &  5.960e+01,5.313e+01,4.293e+01,3.559e+01,3.345e+01,3.517e+01,
     &  3.938e+01,4.437e+01,4.774e+01,4.693e+01,4.244e+01,3.788e+01,
     &  3.560e+01,3.646e+01,4.105e+01,4.730e+01,5.110e+01,5.002e+01,
     &  4.533e+01,4.128e+01,4.008e+01,4.000e+01,3.910e+01,3.733e+01,
     &  3.663e+01,3.779e+01,3.919e+01,4.106e+01,4.521e+01,5.074e+01/
      data (sun2b(i),i=1141,1200)/
     &  5.428e+01,5.501e+01,5.614e+01,6.093e+01,6.738e+01,6.953e+01,
     &  6.516e+01,5.637e+01,4.733e+01,4.432e+01,4.994e+01,5.970e+01,
     &  6.793e+01,7.133e+01,6.990e+01,6.562e+01,5.930e+01,5.402e+01,
     &  5.588e+01,6.515e+01,7.407e+01,7.622e+01,7.351e+01,7.141e+01,
     &  7.085e+01,6.975e+01,6.927e+01,7.138e+01,7.272e+01,6.893e+01,
     &  6.167e+01,5.490e+01,4.742e+01,3.833e+01,3.222e+01,3.124e+01,
     &  3.331e+01,3.536e+01,3.562e+01,3.684e+01,4.155e+01,4.750e+01,
     &  5.118e+01,5.034e+01,4.536e+01,3.834e+01,3.313e+01,3.380e+01,
     &  4.014e+01,4.912e+01,5.538e+01,5.517e+01,5.045e+01,4.651e+01,
     &  4.749e+01,5.188e+01,5.635e+01,5.960e+01,6.158e+01,6.322e+01/
      data (sun2b(i),i=1201,1260)/
     &  6.460e+01,6.410e+01,5.903e+01,5.096e+01,4.763e+01,5.254e+01,
     &  5.888e+01,5.983e+01,5.762e+01,5.673e+01,5.737e+01,5.790e+01,
     &  5.718e+01,5.513e+01,5.295e+01,5.202e+01,5.219e+01,5.204e+01,
     &  5.027e+01,4.659e+01,4.251e+01,4.076e+01,4.189e+01,4.412e+01,
     &  4.654e+01,4.886e+01,5.049e+01,5.192e+01,5.408e+01,5.471e+01,
     &  5.193e+01,4.945e+01,4.987e+01,5.093e+01,5.050e+01,4.862e+01,
     &  4.672e+01,4.607e+01,4.626e+01,4.673e+01,4.801e+01,5.019e+01,
     &  5.242e+01,5.354e+01,5.251e+01,5.138e+01,5.321e+01,5.699e+01,
     &  6.061e+01,6.314e+01,6.400e+01,6.387e+01,6.510e+01,6.939e+01,
     &  7.474e+01,7.818e+01,7.810e+01,7.411e+01,6.737e+01,6.085e+01/
      data (sun2b(i),i=1261,1320)/
     &  5.892e+01,6.268e+01,6.803e+01,6.912e+01,6.460e+01,5.911e+01,
     &  5.600e+01,5.684e+01,6.178e+01,6.587e+01,6.508e+01,6.304e+01,
     &  6.481e+01,6.991e+01,7.484e+01,7.644e+01,7.359e+01,6.885e+01,
     &  6.750e+01,7.268e+01,8.060e+01,8.342e+01,7.896e+01,7.223e+01,
     &  6.674e+01,6.284e+01,6.153e+01,6.357e+01,6.925e+01,7.658e+01,
     &  7.992e+01,7.775e+01,7.394e+01,7.052e+01,6.800e+01,6.634e+01,
     &  6.398e+01,6.110e+01,5.942e+01,5.810e+01,5.574e+01,5.255e+01,
     &  4.808e+01,4.258e+01,3.837e+01,3.730e+01,3.746e+01,3.486e+01,
     &  3.048e+01,2.963e+01,3.473e+01,4.246e+01,4.707e+01,4.585e+01,
     &  4.016e+01,3.429e+01,3.158e+01,3.065e+01,2.905e+01,2.779e+01/
      data (sun2b(i),i=1321,1380)/
     &  3.043e+01,3.757e+01,4.420e+01,4.688e+01,4.785e+01,4.917e+01,
     &  4.918e+01,4.500e+01,3.813e+01,3.506e+01,3.810e+01,4.175e+01,
     &  4.090e+01,3.572e+01,2.888e+01,2.483e+01,2.867e+01,3.965e+01,
     &  5.031e+01,5.572e+01,5.740e+01,5.811e+01,5.941e+01,5.936e+01,
     &  5.342e+01,4.300e+01,3.479e+01,3.370e+01,3.968e+01,4.755e+01,
     &  5.260e+01,5.363e+01,5.100e+01,4.527e+01,3.784e+01,3.103e+01,
     &  2.594e+01,2.280e+01,2.188e+01,2.348e+01,2.786e+01,3.345e+01,
     &  3.732e+01,3.919e+01,4.283e+01,5.040e+01,5.875e+01,6.330e+01,
     &  6.109e+01,5.353e+01,4.605e+01,4.112e+01,3.765e+01,3.630e+01,
     &  4.043e+01,5.089e+01,6.155e+01,6.539e+01,6.268e+01,5.809e+01/
      data (sun2b(i),i=1381,1440)/
     &  5.462e+01,5.133e+01,4.687e+01,4.287e+01,4.055e+01,3.976e+01,
     &  4.022e+01,4.036e+01,3.956e+01,4.067e+01,4.626e+01,5.341e+01,
     &  5.604e+01,5.257e+01,4.667e+01,4.107e+01,3.551e+01,3.123e+01,
     &  3.108e+01,3.596e+01,4.520e+01,5.546e+01,6.180e+01,6.351e+01,
     &  6.185e+01,5.641e+01,4.939e+01,4.637e+01,5.006e+01,5.669e+01,
     &  6.088e+01,6.103e+01,5.811e+01,5.430e+01,5.194e+01,5.051e+01,
     &  4.675e+01,3.915e+01,3.153e+01,2.896e+01,3.097e+01,3.267e+01,
     &  3.157e+01,2.934e+01,2.727e+01,2.518e+01,2.426e+01,2.707e+01,
     &  3.430e+01,4.247e+01,4.723e+01,4.742e+01,4.444e+01,4.054e+01,
     &  3.687e+01,3.302e+01,2.941e+01,2.875e+01,3.408e+01,4.425e+01/
      data (sun2b(i),i=1441,1500)/
     &  5.378e+01,5.797e+01,5.638e+01,5.120e+01,4.531e+01,4.027e+01,
     &  3.590e+01,3.334e+01,3.401e+01,3.686e+01,4.128e+01,4.737e+01,
     &  5.309e+01,5.620e+01,5.563e+01,5.084e+01,4.400e+01,3.877e+01,
     &  3.625e+01,3.638e+01,4.076e+01,5.070e+01,6.337e+01,7.343e+01,
     &  7.642e+01,7.037e+01,5.874e+01,4.703e+01,3.860e+01,3.466e+01,
     &  3.579e+01,4.208e+01,4.997e+01,5.434e+01,5.396e+01,5.229e+01,
     &  5.278e+01,5.557e+01,5.903e+01,6.027e+01,5.625e+01,4.736e+01,
     &  3.806e+01,3.289e+01,3.174e+01,3.173e+01,3.248e+01,3.506e+01,
     &  3.909e+01,4.340e+01,4.813e+01,5.357e+01,5.875e+01,6.360e+01,
     &  6.897e+01,7.342e+01,7.386e+01,6.900e+01,6.056e+01,5.187e+01/
      data (sun2b(i),i=1501,1560)/
     &  4.488e+01,4.206e+01,4.480e+01,4.795e+01,4.688e+01,4.297e+01,
     &  3.929e+01,3.771e+01,3.714e+01,3.522e+01,3.224e+01,3.049e+01,
     &  3.260e+01,4.043e+01,5.110e+01,5.771e+01,5.715e+01,5.299e+01,
     &  5.028e+01,4.999e+01,4.978e+01,4.837e+01,4.642e+01,4.460e+01,
     &  4.273e+01,4.124e+01,4.156e+01,4.381e+01,4.701e+01,4.899e+01,
     &  4.643e+01,4.060e+01,3.784e+01,4.235e+01,5.225e+01,6.053e+01,
     &  6.157e+01,5.680e+01,5.204e+01,5.226e+01,5.708e+01,6.102e+01,
     &  6.071e+01,5.705e+01,5.148e+01,4.635e+01,4.437e+01,4.495e+01,
     &  4.548e+01,4.494e+01,4.383e+01,4.210e+01,3.947e+01,3.683e+01,
     &  3.591e+01,3.636e+01,3.566e+01,3.395e+01,3.369e+01,3.443e+01/
      data (sun2b(i),i=1561,1620)/
     &  3.400e+01,3.265e+01,3.141e+01,3.028e+01,2.941e+01,2.913e+01,
     &  2.933e+01,2.987e+01,3.060e+01,3.131e+01,3.210e+01,3.278e+01,
     &  3.276e+01,3.210e+01,3.198e+01,3.348e+01,3.605e+01,3.917e+01,
     &  4.337e+01,4.724e+01,4.821e+01,4.579e+01,4.159e+01,3.878e+01,
     &  4.075e+01,4.675e+01,5.168e+01,5.260e+01,5.145e+01,5.068e+01,
     &  4.945e+01,4.675e+01,4.537e+01,4.769e+01,5.024e+01,4.896e+01,
     &  4.669e+01,4.860e+01,5.369e+01,5.647e+01,5.434e+01,5.072e+01,
     &  4.988e+01,5.125e+01,5.209e+01,5.276e+01,5.625e+01,6.333e+01,
     &  6.974e+01,7.107e+01,6.835e+01,6.512e+01,6.255e+01,5.919e+01,
     &  5.371e+01,4.816e+01,4.624e+01,4.771e+01,4.966e+01,5.080e+01/
      data (sun2b(i),i=1621,1680)/
     &  5.184e+01,5.454e+01,5.965e+01,6.471e+01,6.584e+01,6.163e+01,
     &  5.557e+01,5.408e+01,5.878e+01,6.489e+01,6.978e+01,7.401e+01,
     &  7.649e+01,7.623e+01,7.475e+01,7.494e+01,7.780e+01,7.962e+01,
     &  7.619e+01,6.719e+01,5.523e+01,4.581e+01,4.314e+01,4.565e+01,
     &  4.947e+01,5.223e+01,5.222e+01,4.889e+01,4.472e+01,4.261e+01,
     &  4.338e+01,4.597e+01,4.812e+01,4.900e+01,4.988e+01,5.071e+01,
     &  4.989e+01,4.832e+01,4.820e+01,5.028e+01,5.383e+01,5.591e+01,
     &  5.482e+01,5.294e+01,5.194e+01,4.944e+01,4.296e+01,3.461e+01,
     &  2.810e+01,2.450e+01,2.420e+01,2.784e+01,3.460e+01,4.162e+01,
     &  4.532e+01,4.544e+01,4.553e+01,4.718e+01,4.576e+01,3.686e+01/
      data (sun2b(i),i=1681,1740)/
     &  2.604e+01,2.057e+01,2.033e+01,2.426e+01,3.086e+01,3.594e+01,
     &  3.671e+01,3.569e+01,3.726e+01,4.086e+01,4.442e+01,4.880e+01,
     &  5.418e+01,5.765e+01,5.843e+01,5.997e+01,6.394e+01,6.682e+01,
     &  6.546e+01,5.948e+01,4.940e+01,3.942e+01,3.418e+01,3.539e+01,
     &  4.288e+01,5.203e+01,5.760e+01,5.909e+01,5.727e+01,5.217e+01,
     &  4.549e+01,3.942e+01,3.558e+01,3.590e+01,4.035e+01,4.673e+01,
     &  5.331e+01,5.878e+01,6.178e+01,5.926e+01,5.003e+01,4.157e+01,
     &  4.052e+01,4.358e+01,4.488e+01,4.275e+01,4.008e+01,3.994e+01,
     &  4.098e+01,3.957e+01,3.496e+01,3.042e+01,3.104e+01,3.869e+01,
     &  4.748e+01,4.983e+01,4.679e+01,4.483e+01,4.655e+01,5.042e+01/
      data (sun2b(i),i=1741,1800)/
     &  5.460e+01,5.766e+01,5.846e+01,5.728e+01,5.562e+01,5.451e+01,
     &  5.334e+01,5.003e+01,4.282e+01,3.364e+01,2.713e+01,2.552e+01,
     &  2.790e+01,3.139e+01,3.213e+01,2.946e+01,2.658e+01,2.596e+01,
     &  2.774e+01,3.117e+01,3.496e+01,3.767e+01,3.864e+01,3.796e+01,
     &  3.673e+01,3.568e+01,3.388e+01,3.085e+01,2.806e+01,2.761e+01,
     &  2.932e+01,2.938e+01,2.539e+01,2.066e+01,1.948e+01,2.230e+01,
     &  2.728e+01,3.247e+01,3.591e+01,3.714e+01,3.790e+01,3.913e+01,
     &  3.978e+01,3.987e+01,4.078e+01,4.232e+01,4.293e+01,4.043e+01,
     &  3.423e+01,2.770e+01,2.388e+01,2.217e+01,2.164e+01,2.259e+01,
     &  2.518e+01,2.902e+01,3.298e+01,3.611e+01,3.858e+01,4.124e+01/
      data (sun2b(i),i=1801,1860)/
     &  4.443e+01,4.694e+01,4.701e+01,4.417e+01,3.966e+01,3.556e+01,
     &  3.284e+01,3.155e+01,3.268e+01,3.696e+01,4.233e+01,4.493e+01,
     &  4.370e+01,4.094e+01,3.797e+01,3.520e+01,3.357e+01,3.334e+01,
     &  3.419e+01,3.635e+01,3.996e+01,4.396e+01,4.716e+01,4.899e+01,
     &  4.898e+01,4.795e+01,4.800e+01,4.989e+01,5.106e+01,4.783e+01,
     &  4.049e+01,3.267e+01,2.680e+01,2.446e+01,2.666e+01,3.193e+01,
     &  3.763e+01,4.135e+01,4.096e+01,3.683e+01,3.211e+01,2.861e+01,
     &  2.648e+01,2.660e+01,2.883e+01,3.088e+01,3.098e+01,3.006e+01,
     &  2.989e+01,3.031e+01,2.997e+01,2.826e+01,2.652e+01,2.707e+01,
     &  3.040e+01,3.454e+01,3.710e+01,3.760e+01,3.725e+01,3.706e+01/
      data (sun2b(i),i=1861,1920)/
     &  3.650e+01,3.417e+01,2.981e+01,2.419e+01,1.851e+01,1.509e+01,
     &  1.504e+01,1.716e+01,2.081e+01,2.568e+01,3.035e+01,3.420e+01,
     &  3.790e+01,4.253e+01,4.783e+01,5.051e+01,4.802e+01,4.262e+01,
     &  3.832e+01,3.737e+01,4.017e+01,4.440e+01,4.613e+01,4.391e+01,
     &  3.840e+01,3.138e+01,2.627e+01,2.508e+01,2.665e+01,2.896e+01,
     &  3.117e+01,3.417e+01,3.805e+01,4.023e+01,3.835e+01,3.274e+01,
     &  2.620e+01,2.186e+01,2.025e+01,2.018e+01,2.173e+01,2.556e+01,
     &  3.032e+01,3.343e+01,3.445e+01,3.496e+01,3.637e+01,3.787e+01,
     &  3.665e+01,3.197e+01,2.592e+01,2.126e+01,2.066e+01,2.466e+01,
     &  3.026e+01,3.402e+01,3.434e+01,3.136e+01,2.693e+01,2.311e+01/
      data (sun2b(i),i=1921,1980)/
     &  2.087e+01,2.068e+01,2.242e+01,2.488e+01,2.678e+01,2.733e+01,
     &  2.654e+01,2.521e+01,2.401e+01,2.294e+01,2.180e+01,2.045e+01,
     &  1.929e+01,1.953e+01,2.174e+01,2.413e+01,2.399e+01,2.156e+01,
     &  1.956e+01,1.892e+01,1.834e+01,1.733e+01,1.655e+01,1.648e+01,
     &  1.721e+01,1.844e+01,1.929e+01,1.898e+01,1.791e+01,1.716e+01,
     &  1.710e+01,1.726e+01,1.692e+01,1.584e+01,1.448e+01,1.368e+01,
     &  1.365e+01,1.391e+01,1.401e+01,1.377e+01,1.346e+01,1.340e+01,
     &  1.355e+01,1.376e+01,1.406e+01,1.443e+01,1.464e+01,1.444e+01,
     &  1.399e+01,1.375e+01,1.393e+01,1.439e+01,1.476e+01,1.482e+01,
     &  1.468e+01,1.445e+01,1.419e+01,1.394e+01,1.375e+01,1.366e+01/
      data (sun2b(i),i=1981,2040)/
     &  1.354e+01,1.331e+01,1.305e+01,1.284e+01,1.270e+01,1.274e+01,
     &  1.281e+01,1.266e+01,1.235e+01,1.210e+01,1.200e+01,1.201e+01,
     &  1.207e+01,1.222e+01,1.244e+01,1.247e+01,1.216e+01,1.173e+01,
     &  1.152e+01,1.162e+01,1.187e+01,1.203e+01,1.195e+01,1.172e+01,
     &  1.140e+01,1.093e+01,1.047e+01,1.020e+01,1.022e+01,1.069e+01,
     &  1.147e+01,1.201e+01,1.208e+01,1.190e+01,1.150e+01,1.089e+01,
     &  1.057e+01,1.085e+01,1.135e+01,1.161e+01,1.141e+01,1.088e+01,
     &  1.038e+01,1.031e+01,1.069e+01,1.124e+01,1.164e+01,1.183e+01,
     &  1.192e+01,1.186e+01,1.167e+01,1.151e+01,1.141e+01,1.130e+01,
     &  1.122e+01,1.114e+01,1.104e+01,1.098e+01,1.100e+01,1.090e+01/
      data (sun2b(i),i=2041,2100)/
     &  1.065e+01,1.056e+01,1.078e+01,1.119e+01,1.160e+01,1.181e+01,
     &  1.158e+01,1.106e+01,1.057e+01,1.034e+01,1.041e+01,1.073e+01,
     &  1.116e+01,1.154e+01,1.165e+01,1.137e+01,1.093e+01,1.052e+01,
     &  9.997e+00,9.378e+00,8.988e+00,9.016e+00,9.413e+00,9.918e+00,
     &  1.028e+01,1.047e+01,1.055e+01,1.057e+01,1.052e+01,1.022e+01,
     &  9.682e+00,9.150e+00,8.706e+00,8.389e+00,8.344e+00,8.624e+00,
     &  9.137e+00,9.718e+00,1.019e+01,1.044e+01,1.046e+01,1.029e+01,
     &  9.977e+00,9.583e+00,9.310e+00,9.319e+00,9.469e+00,9.518e+00,
     &  9.433e+00,9.248e+00,8.820e+00,7.981e+00,7.000e+00,6.486e+00,
     &  6.746e+00,7.541e+00,8.281e+00,8.726e+00,9.068e+00,9.167e+00/
      data (sun2b(i),i=2101,2160)/
     &  8.864e+00,8.464e+00,8.239e+00,8.166e+00,8.188e+00,8.358e+00,
     &  8.649e+00,8.898e+00,9.012e+00,9.066e+00,9.196e+00,9.421e+00,
     &  9.553e+00,9.424e+00,9.129e+00,8.844e+00,8.614e+00,8.424e+00,
     &  8.298e+00,8.260e+00,8.286e+00,8.347e+00,8.453e+00,8.628e+00,
     &  8.831e+00,8.887e+00,8.675e+00,8.331e+00,8.009e+00,7.730e+00,
     &  7.624e+00,7.869e+00,8.273e+00,8.409e+00,8.252e+00,8.091e+00,
     &  8.000e+00,7.937e+00,7.927e+00,7.958e+00,8.049e+00,8.238e+00,
     &  8.426e+00,8.481e+00,8.331e+00,8.026e+00,7.763e+00,7.699e+00,
     &  7.812e+00,7.939e+00,8.018e+00,8.082e+00,8.043e+00,7.892e+00,
     &  7.696e+00,7.497e+00,7.401e+00,7.429e+00,7.449e+00,7.406e+00/
      data (sun2b(i),i=2161,2220)/
     &  7.420e+00,7.526e+00,7.525e+00,7.324e+00,7.126e+00,7.142e+00,
     &  7.334e+00,7.505e+00,7.548e+00,7.532e+00,7.516e+00,7.499e+00,
     &  7.573e+00,7.810e+00,8.059e+00,8.095e+00,7.914e+00,7.698e+00,
     &  7.512e+00,7.214e+00,6.806e+00,6.543e+00,6.523e+00,6.587e+00,
     &  6.561e+00,6.398e+00,6.189e+00,6.059e+00,6.068e+00,6.199e+00,
     &  6.314e+00,6.253e+00,6.093e+00,6.028e+00,6.094e+00,6.303e+00,
     &  6.659e+00,6.940e+00,6.957e+00,6.831e+00,6.737e+00,6.681e+00,
     &  6.656e+00,6.834e+00,7.202e+00,7.401e+00,7.295e+00,7.049e+00,
     &  6.797e+00,6.629e+00,6.587e+00,6.598e+00,6.601e+00,6.650e+00,
     &  6.763e+00,6.785e+00,6.632e+00,6.486e+00,6.520e+00,6.649e+00/
      data (sun2b(i),i=2221,2280)/
     &  6.727e+00,6.723e+00,6.670e+00,6.619e+00,6.598e+00,6.619e+00,
     &  6.711e+00,6.834e+00,6.875e+00,6.825e+00,6.788e+00,6.827e+00,
     &  6.856e+00,6.807e+00,6.838e+00,7.047e+00,7.278e+00,7.414e+00,
     &  7.471e+00,7.462e+00,7.407e+00,7.359e+00,7.337e+00,7.322e+00,
     &  7.294e+00,7.253e+00,7.205e+00,7.134e+00,7.030e+00,6.853e+00,
     &  6.554e+00,6.223e+00,6.014e+00,5.938e+00,5.904e+00,5.857e+00,
     &  5.791e+00,5.733e+00,5.774e+00,5.961e+00,6.187e+00,6.368e+00,
     &  6.500e+00,6.544e+00,6.464e+00,6.271e+00,6.045e+00,5.956e+00,
     &  6.085e+00,6.254e+00,6.245e+00,6.080e+00,5.910e+00,5.845e+00,
     &  5.922e+00,6.106e+00,6.290e+00,6.327e+00,6.110e+00,5.742e+00/
      data (sun2b(i),i=2281,2340)/
     &  5.445e+00,5.298e+00,5.326e+00,5.494e+00,5.682e+00,5.825e+00,
     &  5.893e+00,5.963e+00,6.170e+00,6.452e+00,6.633e+00,6.696e+00,
     &  6.719e+00,6.624e+00,6.311e+00,5.924e+00,5.699e+00,5.665e+00,
     &  5.743e+00,5.879e+00,5.971e+00,5.962e+00,5.967e+00,6.075e+00,
     &  6.254e+00,6.430e+00,6.497e+00,6.456e+00,6.408e+00,6.302e+00,
     &  6.014e+00,5.643e+00,5.396e+00,5.299e+00,5.264e+00,5.223e+00,
     &  5.128e+00,4.931e+00,4.635e+00,4.317e+00,4.015e+00,3.662e+00,
     &  3.291e+00,3.103e+00,3.135e+00,3.199e+00,3.260e+00,3.414e+00,
     &  3.595e+00,3.653e+00,3.630e+00,3.628e+00,3.588e+00,3.429e+00,
     &  3.266e+00,3.228e+00,3.258e+00,3.300e+00,3.454e+00,3.737e+00/
      data (sun2b(i),i=2341,2400)/
     &  4.057e+00,4.356e+00,4.600e+00,4.778e+00,4.885e+00,4.900e+00,
     &  4.839e+00,4.762e+00,4.706e+00,4.698e+00,4.767e+00,4.845e+00,
     &  4.824e+00,4.729e+00,4.686e+00,4.713e+00,4.773e+00,4.871e+00,
     &  4.960e+00,4.953e+00,4.896e+00,4.925e+00,5.074e+00,5.223e+00,
     &  5.251e+00,5.154e+00,5.016e+00,4.888e+00,4.769e+00,4.655e+00,
     &  4.553e+00,4.483e+00,4.466e+00,4.504e+00,4.590e+00,4.703e+00,
     &  4.786e+00,4.833e+00,4.928e+00,5.038e+00,5.006e+00,4.847e+00,
     &  4.683e+00,4.559e+00,4.481e+00,4.431e+00,4.390e+00,4.383e+00,
     &  4.407e+00,4.390e+00,4.297e+00,4.198e+00,4.146e+00,4.108e+00,
     &  4.149e+00,4.390e+00,4.686e+00,4.821e+00,4.794e+00,4.678e+00/
      data (sun2b(i),i=2401,2460)/
     &  4.522e+00,4.396e+00,4.336e+00,4.295e+00,4.226e+00,4.145e+00,
     &  4.106e+00,4.146e+00,4.215e+00,4.255e+00,4.306e+00,4.374e+00,
     &  4.374e+00,4.272e+00,4.139e+00,4.041e+00,3.946e+00,3.813e+00,
     &  3.710e+00,3.734e+00,3.859e+00,3.960e+00,3.953e+00,3.838e+00,
     &  3.671e+00,3.536e+00,3.493e+00,3.537e+00,3.638e+00,3.789e+00,
     &  3.967e+00,4.100e+00,4.142e+00,4.130e+00,4.123e+00,4.162e+00,
     &  4.251e+00,4.295e+00,4.223e+00,4.099e+00,3.968e+00,3.947e+00,
     &  4.196e+00,4.513e+00,4.610e+00,4.513e+00,4.347e+00,4.175e+00,
     &  4.047e+00,3.996e+00,4.008e+00,4.037e+00,4.025e+00,3.938e+00,
     &  3.819e+00,3.753e+00,3.767e+00,3.826e+00,3.865e+00,3.852e+00/
      data (sun2b(i),i=2461,2520)/
     &  3.815e+00,3.803e+00,3.810e+00,3.793e+00,3.766e+00,3.794e+00,
     &  3.869e+00,3.898e+00,3.886e+00,3.892e+00,3.857e+00,3.694e+00,
     &  3.469e+00,3.322e+00,3.282e+00,3.289e+00,3.304e+00,3.322e+00,
     &  3.331e+00,3.333e+00,3.348e+00,3.381e+00,3.419e+00,3.460e+00,
     &  3.497e+00,3.515e+00,3.516e+00,3.516e+00,3.526e+00,3.541e+00,
     &  3.548e+00,3.534e+00,3.534e+00,3.582e+00,3.640e+00,3.645e+00,
     &  3.587e+00,3.511e+00,3.480e+00,3.526e+00,3.624e+00,3.729e+00,
     &  3.802e+00,3.802e+00,3.727e+00,3.658e+00,3.644e+00,3.642e+00,
     &  3.637e+00,3.640e+00,3.641e+00,3.649e+00,3.661e+00,3.625e+00,
     &  3.521e+00,3.402e+00,3.284e+00,3.123e+00,2.948e+00,2.852e+00/
      data (sun2b(i),i=2521,2580)/
     &  2.843e+00,2.864e+00,2.897e+00,2.951e+00,3.003e+00,3.055e+00,
     &  3.111e+00,3.118e+00,3.060e+00,2.985e+00,2.932e+00,2.890e+00,
     &  2.859e+00,2.850e+00,2.845e+00,2.812e+00,2.763e+00,2.742e+00,
     &  2.767e+00,2.802e+00,2.816e+00,2.811e+00,2.813e+00,2.857e+00,
     &  2.966e+00,3.106e+00,3.199e+00,3.213e+00,3.209e+00,3.239e+00,
     &  3.266e+00,3.236e+00,3.117e+00,2.909e+00,2.695e+00,2.532e+00,
     &  2.396e+00,2.295e+00,2.251e+00,2.224e+00,2.181e+00,2.130e+00,
     &  2.148e+00,2.326e+00,2.586e+00,2.723e+00,2.649e+00,2.451e+00,
     &  2.244e+00,2.085e+00,1.989e+00,1.984e+00,2.082e+00,2.223e+00,
     &  2.325e+00,2.355e+00,2.348e+00,2.361e+00,2.430e+00,2.536e+00/
      data (sun2b(i),i=2581,2640)/
     &  2.613e+00,2.622e+00,2.607e+00,2.620e+00,2.646e+00,2.643e+00,
     &  2.582e+00,2.468e+00,2.332e+00,2.240e+00,2.264e+00,2.397e+00,
     &  2.552e+00,2.689e+00,2.851e+00,3.010e+00,3.068e+00,3.006e+00,
     &  2.911e+00,2.861e+00,2.852e+00,2.837e+00,2.789e+00,2.729e+00,
     &  2.714e+00,2.773e+00,2.871e+00,2.954e+00,2.995e+00,2.991e+00,
     &  2.940e+00,2.855e+00,2.773e+00,2.730e+00,2.737e+00,2.765e+00,
     &  2.770e+00,2.737e+00,2.683e+00,2.622e+00,2.566e+00,2.534e+00,
     &  2.535e+00,2.568e+00,2.612e+00,2.631e+00,2.602e+00,2.543e+00,
     &  2.484e+00,2.455e+00,2.468e+00,2.510e+00,2.545e+00,2.553e+00,
     &  2.556e+00,2.589e+00,2.662e+00,2.755e+00,2.847e+00,2.929e+00/
      data (sun2b(i),i=2641,2700)/
     &  2.996e+00,3.017e+00,2.977e+00,2.899e+00,2.800e+00,2.693e+00,
     &  2.617e+00,2.593e+00,2.581e+00,2.528e+00,2.428e+00,2.336e+00,
     &  2.316e+00,2.385e+00,2.487e+00,2.546e+00,2.577e+00,2.681e+00,
     &  2.972e+00,3.576e+00,4.401e+00,5.004e+00,5.030e+00,4.514e+00,
     &  3.707e+00,2.906e+00,2.360e+00,2.142e+00,2.112e+00,2.087e+00,
     &  2.030e+00,2.004e+00,2.036e+00,2.096e+00,2.152e+00,2.188e+00,
     &  2.195e+00,2.186e+00,2.190e+00,2.217e+00,2.257e+00,2.289e+00,
     &  2.278e+00,2.215e+00,2.164e+00,2.234e+00,2.473e+00,2.812e+00,
     &  3.129e+00,3.298e+00,3.221e+00,2.886e+00,2.447e+00,2.144e+00,
     &  2.073e+00,2.139e+00,2.217e+00,2.258e+00,2.265e+00,2.248e+00/
      data (sun2b(i),i=2701,2760)/
     &  2.210e+00,2.166e+00,2.136e+00,2.132e+00,2.144e+00,2.146e+00,
     &  2.140e+00,2.146e+00,2.159e+00,2.158e+00,2.148e+00,2.134e+00,
     &  2.110e+00,2.075e+00,2.040e+00,2.012e+00,1.993e+00,1.980e+00,
     &  1.970e+00,1.961e+00,1.954e+00,1.945e+00,1.931e+00,1.906e+00,
     &  1.873e+00,1.847e+00,1.847e+00,1.869e+00,1.884e+00,1.863e+00,
     &  1.814e+00,1.762e+00,1.719e+00,1.685e+00,1.666e+00,1.671e+00,
     &  1.704e+00,1.752e+00,1.812e+00,1.901e+00,2.012e+00,2.098e+00,
     &  2.138e+00,2.148e+00,2.135e+00,2.109e+00,2.095e+00,2.106e+00,
     &  2.109e+00,2.069e+00,2.000e+00,1.946e+00,1.919e+00,1.909e+00,
     &  1.902e+00,1.898e+00,1.905e+00,1.930e+00,1.959e+00,1.964e+00/
      data (sun2b(i),i=2761,2820)/
     &  1.936e+00,1.902e+00,1.889e+00,1.894e+00,1.902e+00,1.904e+00,
     &  1.898e+00,1.875e+00,1.829e+00,1.772e+00,1.730e+00,1.733e+00,
     &  1.790e+00,1.878e+00,1.955e+00,1.991e+00,1.981e+00,1.943e+00,
     &  1.917e+00,1.922e+00,1.920e+00,1.872e+00,1.783e+00,1.696e+00,
     &  1.639e+00,1.608e+00,1.592e+00,1.582e+00,1.558e+00,1.514e+00,
     &  1.462e+00,1.424e+00,1.425e+00,1.483e+00,1.597e+00,1.741e+00,
     &  1.877e+00,1.978e+00,2.045e+00,2.087e+00,2.091e+00,2.038e+00,
     &  1.957e+00,1.900e+00,1.882e+00,1.866e+00,1.819e+00,1.754e+00,
     &  1.687e+00,1.622e+00,1.573e+00,1.545e+00,1.529e+00,1.531e+00,
     &  1.570e+00,1.641e+00,1.719e+00,1.768e+00,1.755e+00,1.662e+00/
      data (sun2b(i),i=2821,2880)/
     &  1.512e+00,1.359e+00,1.260e+00,1.235e+00,1.261e+00,1.309e+00,
     &  1.359e+00,1.388e+00,1.380e+00,1.348e+00,1.322e+00,1.316e+00,
     &  1.344e+00,1.414e+00,1.495e+00,1.544e+00,1.552e+00,1.528e+00,
     &  1.490e+00,1.461e+00,1.446e+00,1.440e+00,1.440e+00,1.454e+00,
     &  1.476e+00,1.478e+00,1.451e+00,1.423e+00,1.418e+00,1.422e+00,
     &  1.412e+00,1.391e+00,1.378e+00,1.381e+00,1.393e+00,1.409e+00,
     &  1.427e+00,1.442e+00,1.441e+00,1.419e+00,1.386e+00,1.355e+00,
     &  1.326e+00,1.298e+00,1.277e+00,1.273e+00,1.286e+00,1.300e+00,
     &  1.306e+00,1.299e+00,1.282e+00,1.259e+00,1.229e+00,1.187e+00,
     &  1.143e+00,1.118e+00,1.114e+00,1.112e+00,1.101e+00,1.080e+00/
      data (sun2b(i),i=2881,2910)/
     &  1.052e+00,1.028e+00,1.025e+00,1.058e+00,1.120e+00,1.179e+00,
     &  1.206e+00,1.201e+00,1.194e+00,1.200e+00,1.204e+00,1.182e+00,
     &  1.149e+00,1.140e+00,1.154e+00,1.163e+00,1.162e+00,1.159e+00,
     &  1.156e+00,1.157e+00,1.156e+00,1.140e+00,1.104e+00,1.063e+00,
     &  1.034e+00,1.022e+00,1.020e+00,1.019e+00,1.017e+00,1.016e+00/

c sun2a   lowtran7        0-28780 cm-1      20 cm-1           1440
c sun2b   lowtran7    28400-57490 cm-1      10 cm-1           2910

      wnmin=28400.
      wnmax=57490.
      nn1=2910
      do i=1,nn1
        sun(i)=sun2b(2910-i+1)
        wn=wnmin+(wnmax-wnmin)*real(nn1-i)/(nn1-1)
        wlsun(i)=10000./wn
      enddo
      wnmin=0.
      wnmax=28780.
      nn2=1440
      do i=1,nn2
        sun(i+nn1)=sun2a(nn2-i+1)
        wn=wnmin+(wnmax-wnmin)*real(nn2-i)/(nn2-1)
        wlsun(i+nn1)=10000./max(wn,one)
      enddo
      nns=nn1+nn2
      return
      end
c=======================================================================
      subroutine sunmod(wlsun,sun,nns)
      use params, only: kr
      implicit none
      integer :: i, nns
      real(kr) :: wn, wnmax, wnmin, wlsun(*), sun(*), sun3(2494)

c block    source           range           resolution      n_elements
c -----    ------           -----           ----------      ----------
c sun3    modtran3      100-49960 cm-1      20 cm-1           2494

      data (sun3(i),i=1,60)/
     &  2.616e-05,5.360e-05,9.862e-05,1.675e-04,2.676e-04,4.072e-04,
     &  5.959e-04,8.435e-04,1.161e-03,1.561e-03,2.056e-03,2.660e-03,
     &  3.387e-03,4.254e-03,5.277e-03,6.473e-03,7.861e-03,9.461e-03,
     &  1.129e-02,1.338e-02,1.574e-02,1.840e-02,2.138e-02,2.471e-02,
     &  2.840e-02,3.250e-02,3.702e-02,4.200e-02,4.744e-02,5.341e-02,
     &  5.993e-02,6.699e-02,7.469e-02,8.299e-02,9.200e-02,1.017e-01,
     &  1.121e-01,1.234e-01,1.354e-01,1.482e-01,1.620e-01,1.767e-01,
     &  1.924e-01,2.090e-01,2.268e-01,2.455e-01,2.653e-01,2.864e-01,
     &  3.087e-01,3.323e-01,3.571e-01,3.832e-01,4.111e-01,4.403e-01,
     &  4.709e-01,5.029e-01,5.366e-01,5.721e-01,6.095e-01,6.483e-01/
      data (sun3(i),i=61,120)/
     &  6.888e-01,7.313e-01,7.738e-01,8.219e-01,8.701e-01,9.209e-01,
     &  9.728e-01,1.028e+00,1.085e+00,1.143e+00,1.205e+00,1.269e+00,
     &  1.334e+00,1.403e+00,1.473e+00,1.545e+00,1.623e+00,1.697e+00,
     &  1.780e+00,1.861e+00,1.942e+00,2.036e+00,2.125e+00,2.218e+00,
     &  2.311e+00,2.410e+00,2.511e+00,2.613e+00,2.724e+00,2.833e+00,
     &  2.950e+00,3.064e+00,3.188e+00,3.305e+00,3.446e+00,3.570e+00,
     &  3.708e+00,3.854e+00,3.995e+00,4.136e+00,4.285e+00,4.433e+00,
     &  4.562e+00,4.727e+00,4.871e+00,5.046e+00,5.262e+00,5.429e+00,
     &  5.615e+00,5.912e+00,6.135e+00,6.392e+00,6.773e+00,6.996e+00,
     &  7.194e+00,7.452e+00,7.702e+00,7.942e+00,8.132e+00,8.433e+00/
      data (sun3(i),i=121,180)/
     &  8.682e+00,8.981e+00,9.277e+00,9.561e+00,9.832e+00,1.015e+01,
     &  1.044e+01,1.071e+01,1.099e+01,1.126e+01,1.160e+01,1.202e+01,
     &  1.238e+01,1.279e+01,1.313e+01,1.347e+01,1.385e+01,1.425e+01,
     &  1.464e+01,1.502e+01,1.543e+01,1.579e+01,1.616e+01,1.665e+01,
     &  1.709e+01,1.754e+01,1.788e+01,1.826e+01,1.892e+01,1.939e+01,
     &  1.989e+01,2.037e+01,2.090e+01,2.138e+01,2.187e+01,2.239e+01,
     &  2.291e+01,2.350e+01,2.405e+01,2.451e+01,2.505e+01,2.579e+01,
     &  2.634e+01,2.700e+01,2.760e+01,2.817e+01,2.878e+01,2.947e+01,
     &  3.010e+01,3.050e+01,3.142e+01,3.210e+01,3.282e+01,3.350e+01,
     &  3.423e+01,3.481e+01,3.531e+01,3.637e+01,3.704e+01,3.783e+01/
      data (sun3(i),i=181,240)/
     &  3.860e+01,3.937e+01,3.996e+01,4.046e+01,4.133e+01,4.156e+01,
     &  4.284e+01,4.417e+01,4.461e+01,4.582e+01,4.649e+01,4.738e+01,
     &  4.823e+01,4.949e+01,4.999e+01,5.126e+01,5.166e+01,5.270e+01,
     &  5.403e+01,5.429e+01,5.603e+01,5.703e+01,5.750e+01,5.929e+01,
     &  6.004e+01,6.080e+01,6.217e+01,6.299e+01,6.477e+01,6.565e+01,
     &  6.636e+01,6.828e+01,6.955e+01,7.028e+01,7.174e+01,7.342e+01,
     &  7.435e+01,7.547e+01,7.690e+01,7.851e+01,7.960e+01,8.073e+01,
     &  8.188e+01,8.337e+01,8.455e+01,8.548e+01,8.388e+01,8.903e+01,
     &  9.049e+01,9.183e+01,9.305e+01,9.434e+01,9.624e+01,9.754e+01,
     &  9.832e+01,1.002e+02,1.022e+02,1.031e+02,1.050e+02,1.068e+02/
      data (sun3(i),i=241,300)/
     &  1.079e+02,1.091e+02,1.110e+02,1.137e+02,1.151e+02,1.161e+02,
     &  1.179e+02,1.188e+02,1.198e+02,1.230e+02,1.240e+02,1.251e+02,
     &  1.209e+02,1.291e+02,1.299e+02,1.335e+02,1.357e+02,1.372e+02,
     &  1.386e+02,1.373e+02,1.416e+02,1.409e+02,1.381e+02,1.469e+02,
     &  1.496e+02,1.534e+02,1.525e+02,1.566e+02,1.586e+02,1.597e+02,
     &  1.551e+02,1.651e+02,1.670e+02,1.699e+02,1.722e+02,1.742e+02,
     &  1.752e+02,1.789e+02,1.810e+02,1.821e+02,1.851e+02,1.847e+02,
     &  1.874e+02,1.813e+02,1.922e+02,1.946e+02,1.992e+02,2.003e+02,
     &  2.034e+02,2.056e+02,2.103e+02,2.102e+02,2.065e+02,2.106e+02,
     &  2.174e+02,2.200e+02,2.246e+02,2.266e+02,2.268e+02,2.266e+02/
      data (sun3(i),i=301,360)/
     &  2.236e+02,2.355e+02,2.401e+02,2.393e+02,2.363e+02,2.408e+02,
     &  2.434e+02,2.464e+02,2.490e+02,2.479e+02,2.454e+02,2.538e+02,
     &  2.494e+02,2.602e+02,2.603e+02,2.664e+02,2.656e+02,2.683e+02,
     &  2.712e+02,2.753e+02,2.748e+02,2.787e+02,2.791e+02,2.815e+02,
     &  2.844e+02,2.852e+02,2.877e+02,2.798e+02,2.883e+02,2.928e+02,
     &  2.983e+02,2.898e+02,2.998e+02,3.005e+02,3.003e+02,3.058e+02,
     &  3.080e+02,3.095e+02,3.128e+02,3.122e+02,3.177e+02,3.178e+02,
     &  3.161e+02,3.260e+02,3.263e+02,3.245e+02,3.261e+02,3.317e+02,
     &  3.368e+02,3.350e+02,3.419e+02,3.424e+02,3.439e+02,3.496e+02,
     &  3.515e+02,3.534e+02,3.553e+02,3.570e+02,3.580e+02,3.626e+02/
      data (sun3(i),i=361,420)/
     &  3.608e+02,3.661e+02,3.684e+02,3.691e+02,3.716e+02,3.755e+02,
     &  3.792e+02,3.820e+02,3.808e+02,3.857e+02,3.868e+02,3.864e+02,
     &  3.926e+02,3.953e+02,3.973e+02,3.951e+02,3.974e+02,4.050e+02,
     &  4.067e+02,4.071e+02,4.118e+02,4.145e+02,4.156e+02,4.172e+02,
     &  4.172e+02,3.925e+02,4.252e+02,4.295e+02,4.324e+02,4.321e+02,
     &  4.345e+02,4.379e+02,4.412e+02,4.429e+02,4.450e+02,4.501e+02,
     &  4.513e+02,4.513e+02,4.551e+02,4.610e+02,4.613e+02,4.633e+02,
     &  4.643e+02,4.702e+02,4.716e+02,4.748e+02,4.787e+02,4.798e+02,
     &  4.759e+02,4.728e+02,4.782e+02,4.883e+02,4.758e+02,4.890e+02,
     &  4.985e+02,4.974e+02,4.925e+02,5.032e+02,4.967e+02,5.067e+02/
      data (sun3(i),i=421,480)/
     &  5.043e+02,5.155e+02,5.195e+02,5.145e+02,5.168e+02,5.220e+02,
     &  5.154e+02,5.328e+02,5.352e+02,5.348e+02,5.380e+02,5.442e+02,
     &  5.393e+02,5.367e+02,5.405e+02,5.500e+02,5.499e+02,5.546e+02,
     &  5.535e+02,5.538e+02,5.655e+02,5.632e+02,5.667e+02,5.714e+02,
     &  5.733e+02,5.755e+02,5.800e+02,5.817e+02,5.851e+02,5.802e+02,
     &  5.812e+02,5.811e+02,5.676e+02,6.020e+02,6.002e+02,5.986e+02,
     &  6.136e+02,5.933e+02,6.165e+02,6.227e+02,6.190e+02,6.219e+02,
     &  6.277e+02,6.196e+02,6.351e+02,6.417e+02,6.398e+02,6.411e+02,
     &  6.543e+02,6.575e+02,6.563e+02,6.618e+02,6.636e+02,6.595e+02,
     &  6.708e+02,6.734e+02,6.752e+02,6.738e+02,6.779e+02,6.821e+02/
      data (sun3(i),i=481,540)/
     &  6.895e+02,6.921e+02,6.972e+02,6.996e+02,6.967e+02,6.977e+02,
     &  7.032e+02,7.055e+02,7.114e+02,7.164e+02,7.207e+02,7.219e+02,
     &  6.953e+02,7.163e+02,7.332e+02,7.349e+02,7.385e+02,7.406e+02,
     &  7.446e+02,7.453e+02,7.400e+02,7.501e+02,7.525e+02,7.556e+02,
     &  7.620e+02,7.627e+02,7.606e+02,7.633e+02,7.630e+02,7.749e+02,
     &  7.766e+02,7.803e+02,7.783e+02,7.758e+02,7.792e+02,7.886e+02,
     &  7.940e+02,7.932e+02,7.923e+02,7.738e+02,8.072e+02,8.093e+02,
     &  8.181e+02,8.155e+02,8.121e+02,8.103e+02,7.890e+02,8.237e+02,
     &  8.339e+02,8.325e+02,8.386e+02,8.418e+02,8.482e+02,8.440e+02,
     &  8.558e+02,8.305e+02,8.363e+02,8.231e+02,8.442e+02,8.718e+02/
      data (sun3(i),i=541,600)/
     &  8.636e+02,8.772e+02,8.843e+02,8.800e+02,8.759e+02,8.658e+02,
     &  8.900e+02,8.899e+02,9.093e+02,8.799e+02,8.652e+02,9.115e+02,
     &  9.090e+02,9.277e+02,9.107e+02,9.087e+02,9.335e+02,9.345e+02,
     &  9.326e+02,8.978e+02,9.425e+02,9.270e+02,9.424e+02,9.147e+02,
     &  9.525e+02,9.536e+02,9.335e+02,9.519e+02,9.660e+02,9.639e+02,
     &  9.387e+02,9.488e+02,8.045e+02,9.736e+02,9.908e+02,9.743e+02,
     &  9.846e+02,9.786e+02,1.002e+03,9.952e+02,7.876e+02,9.779e+02,
     &  9.915e+02,9.268e+02,1.020e+03,1.008e+03,1.027e+03,1.017e+03,
     &  1.015e+03,1.034e+03,1.039e+03,1.019e+03,1.049e+03,1.047e+03,
     &  1.036e+03,1.025e+03,1.055e+03,1.048e+03,1.065e+03,1.075e+03/
      data (sun3(i),i=601,660)/
     &  1.076e+03,1.072e+03,1.060e+03,1.069e+03,1.062e+03,1.060e+03,
     &  1.072e+03,1.097e+03,1.100e+03,1.109e+03,1.112e+03,1.111e+03,
     &  1.103e+03,1.099e+03,1.114e+03,1.134e+03,1.105e+03,1.137e+03,
     &  1.126e+03,1.145e+03,1.128e+03,1.151e+03,1.153e+03,1.148e+03,
     &  1.125e+03,1.125e+03,1.153e+03,1.159e+03,1.173e+03,1.176e+03,
     &  1.181e+03,1.173e+03,1.187e+03,1.173e+03,1.198e+03,1.192e+03,
     &  1.196e+03,1.208e+03,1.199e+03,1.219e+03,1.179e+03,1.217e+03,
     &  1.200e+03,1.203e+03,1.220e+03,1.226e+03,1.235e+03,1.192e+03,
     &  1.207e+03,1.256e+03,1.257e+03,1.220e+03,1.250e+03,1.256e+03,
     &  1.243e+03,1.253e+03,1.261e+03,1.265e+03,1.284e+03,1.256e+03/
      data (sun3(i),i=661,720)/
     &  1.274e+03,1.273e+03,1.257e+03,1.287e+03,1.289e+03,1.291e+03,
     &  1.273e+03,1.294e+03,1.279e+03,1.260e+03,1.270e+03,1.293e+03,
     &  1.283e+03,1.326e+03,1.312e+03,1.315e+03,1.340e+03,1.308e+03,
     &  1.335e+03,1.312e+03,1.330e+03,1.287e+03,1.347e+03,1.343e+03,
     &  1.353e+03,1.356e+03,1.369e+03,1.332e+03,1.359e+03,1.344e+03,
     &  1.315e+03,1.341e+03,1.382e+03,1.353e+03,1.343e+03,1.393e+03,
     &  1.374e+03,1.384e+03,1.387e+03,1.416e+03,1.396e+03,1.413e+03,
     &  1.393e+03,1.416e+03,1.432e+03,1.422e+03,1.421e+03,1.385e+03,
     &  1.418e+03,1.419e+03,1.454e+03,1.427e+03,1.452e+03,1.450e+03,
     &  1.453e+03,1.453e+03,1.462e+03,1.459e+03,1.455e+03,1.480e+03/
      data (sun3(i),i=721,780)/
     &  1.488e+03,1.469e+03,1.492e+03,1.480e+03,1.449e+03,1.489e+03,
     &  1.479e+03,1.500e+03,1.500e+03,1.490e+03,1.504e+03,1.506e+03,
     &  1.508e+03,1.518e+03,1.513e+03,1.521e+03,1.519e+03,1.530e+03,
     &  1.530e+03,1.519e+03,1.526e+03,1.533e+03,1.554e+03,1.554e+03,
     &  1.552e+03,1.523e+03,1.563e+03,1.547e+03,1.550e+03,1.553e+03,
     &  1.583e+03,1.565e+03,1.586e+03,1.526e+03,1.579e+03,1.543e+03,
     &  1.435e+03,1.230e+03,1.547e+03,1.564e+03,1.595e+03,1.605e+03,
     &  1.590e+03,1.618e+03,1.545e+03,1.536e+03,1.614e+03,1.605e+03,
     &  1.573e+03,1.598e+03,1.616e+03,1.598e+03,1.611e+03,1.634e+03,
     &  1.609e+03,1.605e+03,1.614e+03,1.630e+03,1.644e+03,1.675e+03/
      data (sun3(i),i=781,840)/
     &  1.638e+03,1.606e+03,1.647e+03,1.644e+03,1.618e+03,1.687e+03,
     &  1.584e+03,1.655e+03,1.645e+03,1.683e+03,1.703e+03,1.684e+03,
     &  1.690e+03,1.672e+03,1.587e+03,1.655e+03,1.623e+03,1.648e+03,
     &  1.712e+03,1.658e+03,1.734e+03,1.715e+03,1.684e+03,1.729e+03,
     &  1.712e+03,1.645e+03,1.588e+03,1.678e+03,1.695e+03,1.643e+03,
     &  1.676e+03,1.668e+03,1.735e+03,1.739e+03,1.702e+03,1.756e+03,
     &  1.746e+03,1.738e+03,1.802e+03,1.756e+03,1.783e+03,1.802e+03,
     &  1.806e+03,1.809e+03,1.762e+03,1.743e+03,1.758e+03,1.758e+03,
     &  1.782e+03,1.796e+03,1.786e+03,1.799e+03,1.816e+03,1.827e+03,
     &  1.817e+03,1.782e+03,1.828e+03,1.826e+03,1.777e+03,1.840e+03/
      data (sun3(i),i=841,900)/
     &  1.788e+03,1.820e+03,1.834e+03,1.642e+03,1.648e+03,1.816e+03,
     &  1.864e+03,1.845e+03,1.760e+03,1.838e+03,1.855e+03,1.874e+03,
     &  1.862e+03,1.886e+03,1.874e+03,1.862e+03,1.834e+03,1.878e+03,
     &  1.851e+03,1.802e+03,1.848e+03,1.859e+03,1.875e+03,1.891e+03,
     &  1.803e+03,1.876e+03,1.891e+03,1.911e+03,1.933e+03,1.895e+03,
     &  1.808e+03,1.811e+03,1.842e+03,1.877e+03,1.835e+03,1.861e+03,
     &  1.913e+03,1.844e+03,1.695e+03,1.853e+03,1.898e+03,1.868e+03,
     &  1.869e+03,1.936e+03,1.832e+03,1.772e+03,1.966e+03,1.856e+03,
     &  1.834e+03,1.789e+03,1.792e+03,1.881e+03,1.826e+03,1.832e+03,
     &  1.882e+03,1.873e+03,1.928e+03,1.919e+03,1.884e+03,1.836e+03/
      data (sun3(i),i=901,960)/
     &  1.905e+03,1.944e+03,1.852e+03,1.869e+03,1.874e+03,1.870e+03,
     &  1.886e+03,1.883e+03,1.784e+03,1.945e+03,1.906e+03,1.866e+03,
     &  2.008e+03,1.862e+03,1.984e+03,1.871e+03,1.841e+03,1.908e+03,
     &  1.879e+03,1.839e+03,1.742e+03,1.822e+03,1.897e+03,1.957e+03,
     &  1.913e+03,2.002e+03,1.799e+03,1.910e+03,2.024e+03,1.954e+03,
     &  1.870e+03,1.780e+03,1.853e+03,1.694e+03,1.869e+03,1.932e+03,
     &  2.035e+03,2.026e+03,1.871e+03,1.908e+03,2.022e+03,1.882e+03,
     &  1.964e+03,1.777e+03,1.548e+03,1.731e+03,1.895e+03,1.950e+03,
     &  1.917e+03,2.002e+03,1.849e+03,1.749e+03,1.967e+03,1.864e+03,
     &  1.984e+03,1.736e+03,1.954e+03,1.783e+03,1.773e+03,1.450e+03/
      data (sun3(i),i=961,1020)/
     &  1.883e+03,1.551e+03,1.453e+03,1.813e+03,1.994e+03,1.858e+03,
     &  1.854e+03,1.853e+03,1.855e+03,1.940e+03,1.821e+03,1.956e+03,
     &  2.079e+03,1.982e+03,1.940e+03,1.864e+03,2.056e+03,2.043e+03,
     &  1.902e+03,1.846e+03,1.950e+03,1.976e+03,2.007e+03,2.095e+03,
     &  2.067e+03,1.899e+03,1.970e+03,1.782e+03,1.980e+03,1.921e+03,
     &  1.883e+03,1.896e+03,1.826e+03,1.933e+03,1.817e+03,1.904e+03,
     &  1.971e+03,2.043e+03,1.701e+03,1.918e+03,2.093e+03,1.947e+03,
     &  2.036e+03,1.940e+03,2.122e+03,2.134e+03,2.066e+03,1.939e+03,
     &  1.960e+03,2.112e+03,1.978e+03,1.692e+03,2.074e+03,1.826e+03,
     &  2.060e+03,2.111e+03,2.142e+03,1.782e+03,1.950e+03,2.053e+03/
      data (sun3(i),i=1021,1080)/
     &  1.944e+03,1.821e+03,1.932e+03,1.422e+03,1.792e+03,1.975e+03,
     &  2.007e+03,2.026e+03,1.999e+03,2.095e+03,2.090e+03,2.096e+03,
     &  2.156e+03,2.103e+03,2.064e+03,2.069e+03,2.068e+03,2.139e+03,
     &  2.077e+03,1.972e+03,2.155e+03,2.180e+03,2.155e+03,2.122e+03,
     &  2.006e+03,2.114e+03,2.065e+03,2.182e+03,2.096e+03,2.127e+03,
     &  2.013e+03,2.049e+03,1.991e+03,2.172e+03,2.127e+03,2.039e+03,
     &  2.109e+03,1.982e+03,1.994e+03,2.058e+03,2.160e+03,2.077e+03,
     &  2.084e+03,2.027e+03,2.220e+03,2.086e+03,1.889e+03,2.053e+03,
     &  2.162e+03,2.042e+03,1.953e+03,2.039e+03,2.161e+03,1.941e+03,
     &  2.095e+03,1.999e+03,2.096e+03,1.983e+03,2.085e+03,2.082e+03/
      data (sun3(i),i=1081,1140)/
     &  2.112e+03,2.002e+03,2.077e+03,1.993e+03,2.120e+03,2.017e+03,
     &  1.869e+03,2.139e+03,2.088e+03,2.201e+03,2.035e+03,2.193e+03,
     &  1.990e+03,2.160e+03,1.951e+03,2.020e+03,2.062e+03,1.927e+03,
     &  1.922e+03,1.784e+03,2.011e+03,2.039e+03,2.137e+03,2.166e+03,
     &  2.242e+03,2.191e+03,2.135e+03,1.894e+03,1.964e+03,2.069e+03,
     &  1.913e+03,2.091e+03,2.059e+03,1.914e+03,2.001e+03,1.814e+03,
     &  1.783e+03,1.703e+03,1.915e+03,2.018e+03,1.790e+03,2.047e+03,
     &  1.767e+03,1.838e+03,1.930e+03,1.858e+03,2.124e+03,1.644e+03,
     &  2.019e+03,1.787e+03,1.647e+03,1.700e+03,1.783e+03,1.752e+03,
     &  1.634e+03,1.187e+03,1.864e+03,1.751e+03,1.724e+03,1.737e+03/
      data (sun3(i),i=1141,1200)/
     &  2.002e+03,1.917e+03,1.827e+03,1.669e+03,1.568e+03,1.834e+03,
     &  1.518e+03,1.219e+03,1.926e+03,1.877e+03,1.690e+03,1.182e+03,
     &  1.956e+03,1.933e+03,1.284e+03,8.052e+02,1.086e+03,9.452e+02,
     &  1.145e+03,1.526e+03,1.192e+03,1.498e+03,1.667e+03,1.576e+03,
     &  1.515e+03,1.254e+03,1.602e+03,1.841e+03,1.509e+03,1.805e+03,
     &  1.577e+03,1.566e+03,1.627e+03,1.814e+03,1.623e+03,1.554e+03,
     &  1.743e+03,1.616e+03,1.227e+03,1.923e+03,1.812e+03,1.572e+03,
     &  1.848e+03,1.978e+03,1.688e+03,1.446e+03,1.410e+03,1.780e+03,
     &  1.822e+03,1.558e+03,1.759e+03,1.856e+03,1.663e+03,1.530e+03,
     &  1.777e+03,1.944e+03,1.931e+03,1.905e+03,1.656e+03,1.954e+03/
      data (sun3(i),i=1201,1260)/
     &  1.836e+03,1.583e+03,1.752e+03,1.960e+03,1.672e+03,1.668e+03,
     &  1.729e+03,1.854e+03,1.915e+03,1.920e+03,1.893e+03,1.877e+03,
     &  1.851e+03,1.583e+03,1.310e+03,1.619e+03,1.943e+03,1.964e+03,
     &  1.933e+03,1.559e+03,1.932e+03,1.421e+03,1.840e+03,1.441e+03,
     &  1.633e+03,1.266e+03,1.996e+03,1.645e+03,1.640e+03,1.914e+03,
     &  1.521e+03,1.420e+03,1.697e+03,1.998e+03,1.573e+03,1.522e+03,
     &  1.889e+03,1.730e+03,1.811e+03,1.671e+03,1.828e+03,1.749e+03,
     &  1.840e+03,1.497e+03,1.957e+03,1.658e+03,1.594e+03,1.896e+03,
     &  1.756e+03,1.631e+03,1.444e+03,1.691e+03,1.382e+03,8.553e+02,
     &  4.010e+02,3.110e+02,8.529e+02,1.101e+03,1.349e+03,1.541e+03/
      data (sun3(i),i=1261,1320)/
     &  1.358e+03,1.272e+03,1.161e+03,1.074e+03,8.786e+02,3.822e+02,
     &  2.173e+02,5.268e+02,9.261e+02,1.253e+03,1.066e+03,1.399e+03,
     &  1.506e+03,1.544e+03,1.582e+03,1.102e+03,1.306e+03,1.203e+03,
     &  1.403e+03,1.399e+03,9.520e+02,7.840e+02,1.046e+03,9.084e+02,
     &  7.732e+02,1.194e+03,6.025e+02,9.981e+02,1.144e+03,1.157e+03,
     &  4.782e+02,8.259e+02,1.030e+03,9.441e+02,1.151e+03,1.107e+03,
     &  8.535e+02,5.215e+02,7.978e+02,4.610e+02,8.221e+02,6.668e+02,
     &  7.040e+02,7.287e+02,1.121e+03,8.432e+02,1.049e+03,1.427e+03,
     &  1.115e+03,1.509e+03,1.312e+03,6.955e+02,1.074e+03,1.133e+03,
     &  1.021e+03,1.020e+03,1.725e+03,1.730e+03,1.467e+03,1.406e+03/
      data (sun3(i),i=1321,1380)/
     &  1.436e+03,1.132e+03,9.268e+02,1.102e+03,1.193e+03,7.077e+02,
     &  1.646e+03,1.241e+03,5.566e+02,6.371e+02,7.420e+02,1.420e+03,
     &  9.540e+02,4.338e+02,8.207e+02,1.289e+03,1.063e+03,1.279e+03,
     &  7.421e+02,7.780e+02,1.525e+03,1.708e+03,1.203e+03,8.153e+02,
     &  9.097e+02,1.317e+03,1.352e+03,1.334e+03,1.235e+03,1.418e+03,
     &  1.152e+03,9.400e+02,1.236e+03,1.122e+03,1.159e+03,1.363e+03,
     &  1.387e+03,1.103e+03,1.227e+03,1.348e+03,1.342e+03,1.454e+03,
     &  1.448e+03,1.483e+03,1.161e+03,7.835e+02,9.843e+02,9.535e+02,
     &  1.307e+03,1.182e+03,9.650e+02,8.898e+02,8.146e+02,1.513e+03,
     &  1.224e+03,1.119e+03,5.361e+02,1.182e+03,1.276e+03,1.034e+03/
      data (sun3(i),i=1381,1440)/
     &  6.753e+02,1.069e+03,9.200e+02,1.243e+03,1.355e+03,1.119e+03,
     &  1.181e+03,1.125e+03,7.264e+02,3.851e+02,7.448e+02,5.304e+02,
     &  8.958e+02,9.643e+02,8.306e+02,5.733e+02,9.510e+02,5.896e+02,
     &  1.246e+03,1.098e+03,1.125e+03,1.083e+03,1.063e+03,1.326e+03,
     &  1.075e+03,1.334e+03,1.190e+03,1.331e+03,9.672e+02,1.218e+03,
     &  1.075e+03,1.229e+03,9.998e+02,7.448e+02,9.717e+02,1.197e+03,
     &  1.305e+03,9.595e+02,8.832e+02,1.129e+03,1.162e+03,1.184e+03,
     &  1.342e+03,1.263e+03,9.291e+02,9.784e+02,9.976e+02,8.360e+02,
     &  1.326e+03,1.304e+03,1.144e+03,1.307e+03,1.108e+03,5.392e+02,
     &  1.027e+03,1.248e+03,1.085e+03,8.488e+02,1.214e+03,9.699e+02/
      data (sun3(i),i=1441,1500)/
     &  9.096e+02,1.369e+03,1.118e+03,9.123e+02,1.176e+03,1.046e+03,
     &  8.571e+02,6.842e+02,5.575e+02,9.405e+02,1.046e+03,8.434e+02,
     &  9.864e+02,9.891e+02,7.194e+02,6.664e+02,9.103e+02,1.127e+03,
     &  7.969e+02,6.468e+02,8.385e+02,8.230e+02,8.359e+02,9.173e+02,
     &  8.132e+02,1.177e+03,1.016e+03,9.794e+02,5.923e+02,8.344e+02,
     &  9.039e+02,9.497e+02,8.562e+02,8.964e+02,6.116e+02,1.092e+03,
     &  9.351e+02,6.175e+02,5.049e+02,6.799e+02,8.828e+02,9.557e+02,
     &  4.216e+02,3.581e+02,1.047e+03,1.029e+03,1.087e+03,8.082e+02,
     &  9.585e+02,1.027e+03,1.031e+03,9.652e+02,9.996e+02,7.608e+02,
     &  1.049e+03,1.053e+03,9.285e+02,1.007e+03,1.213e+03,8.617e+02/
      data (sun3(i),i=1501,1560)/
     &  9.966e+02,1.115e+03,1.092e+03,7.323e+02,9.399e+02,9.989e+02,
     &  9.204e+02,7.284e+02,1.037e+03,1.177e+03,1.178e+03,9.469e+02,
     &  1.137e+03,1.201e+03,1.149e+03,9.677e+02,8.670e+02,9.166e+02,
     &  9.096e+02,9.328e+02,9.598e+02,1.088e+03,1.204e+03,9.403e+02,
     &  1.045e+03,1.137e+03,1.043e+03,1.217e+03,9.163e+02,1.141e+03,
     &  1.146e+03,8.977e+02,8.804e+02,8.169e+02,6.701e+02,1.193e+03,
     &  7.830e+02,9.112e+02,5.866e+02,7.196e+02,5.980e+02,7.290e+02,
     &  7.352e+02,5.010e+02,4.252e+02,7.595e+02,6.654e+02,7.978e+02,
     &  8.962e+02,7.806e+02,5.683e+02,6.081e+02,5.947e+02,9.489e+02,
     &  9.340e+02,1.035e+03,7.730e+02,6.878e+02,1.069e+03,6.922e+02/
      data (sun3(i),i=1561,1620)/
     &  6.132e+02,4.768e+02,1.128e+03,8.246e+02,6.237e+02,5.390e+02,
     &  5.126e+02,5.142e+02,1.013e+03,1.022e+03,9.397e+02,1.115e+03,
     &  7.794e+02,6.980e+02,5.552e+02,8.115e+02,4.396e+02,4.824e+02,
     &  4.740e+02,8.546e+02,9.784e+02,8.321e+02,9.878e+02,6.204e+02,
     &  6.790e+02,3.038e+02,9.230e+02,8.131e+02,1.036e+03,8.095e+02,
     &  6.136e+02,7.761e+02,7.880e+02,8.014e+02,7.458e+02,8.269e+02,
     &  6.337e+02,6.414e+02,6.504e+02,8.785e+02,8.504e+02,1.136e+03,
     &  9.529e+02,1.086e+03,9.323e+02,7.545e+02,7.780e+02,2.312e+02,
     &  5.092e+02,5.440e+02,7.925e+02,3.707e+02,5.603e+02,7.485e+02,
     &  5.859e+02,9.702e+02,7.035e+02,5.693e+02,9.062e+02,6.349e+02/
      data (sun3(i),i=1621,1680)/
     &  6.152e+02,8.680e+02,6.633e+02,8.073e+02,6.675e+02,4.866e+02,
     &  5.959e+02,7.526e+02,7.618e+02,3.963e+02,4.099e+02,7.378e+02,
     &  7.068e+02,8.668e+02,7.861e+02,4.051e+02,8.343e+02,7.763e+02,
     &  4.941e+02,6.495e+02,5.904e+02,5.761e+02,7.947e+02,8.966e+02,
     &  6.009e+02,6.533e+02,8.680e+02,3.771e+02,6.289e+02,4.458e+02,
     &  2.255e+02,3.380e+02,5.338e+02,6.043e+02,5.928e+02,5.851e+02,
     &  4.306e+02,3.351e+02,5.880e+02,4.450e+02,3.235e+02,2.487e+02,
     &  5.846e+02,6.711e+02,3.432e+02,5.508e+02,6.835e+02,6.591e+02,
     &  4.658e+02,4.033e+02,2.514e+02,4.175e+02,5.239e+02,6.458e+02,
     &  8.497e+02,6.651e+02,3.887e+02,5.942e+02,3.061e+02,5.717e+02/
      data (sun3(i),i=1681,1740)/
     &  2.725e+02,4.818e+02,7.439e+02,6.226e+02,6.662e+02,5.715e+02,
     &  6.673e+02,3.917e+02,7.972e+02,6.082e+02,3.719e+02,4.311e+02,
     &  7.118e+02,5.424e+02,5.332e+02,6.237e+02,4.479e+02,3.117e+02,
     &  7.760e+02,5.962e+02,7.501e+02,6.620e+02,4.200e+02,4.756e+02,
     &  5.516e+02,5.493e+02,6.381e+02,6.239e+02,5.589e+02,6.256e+02,
     &  6.069e+02,7.299e+02,6.370e+02,6.039e+02,6.705e+02,6.279e+02,
     &  7.775e+02,6.866e+02,7.114e+02,6.659e+02,5.818e+02,6.354e+02,
     &  4.635e+02,4.454e+02,4.498e+02,4.285e+02,4.388e+02,4.641e+02,
     &  4.066e+02,2.264e+02,1.395e+02,2.680e+02,4.347e+02,3.362e+02,
     &  4.051e+02,4.642e+02,5.097e+02,3.298e+02,3.507e+02,4.057e+02/
      data (sun3(i),i=1741,1800)/
     &  4.052e+02,3.658e+02,3.675e+02,3.581e+02,2.417e+02,1.870e+02,
     &  1.356e+02,6.816e+01,5.988e+01,1.213e+02,1.771e+02,2.441e+02,
     &  3.429e+02,2.237e+02,3.777e+02,3.196e+02,3.687e+02,3.892e+02,
     &  2.577e+02,4.301e+02,3.660e+02,3.880e+02,3.655e+02,3.284e+02,
     &  2.890e+02,3.682e+02,3.186e+02,3.097e+02,3.144e+02,2.747e+02,
     &  2.556e+02,2.859e+02,1.762e+02,1.761e+02,1.817e+02,1.524e+02,
     &  1.149e+02,8.834e+01,7.683e+01,7.632e+01,7.916e+01,7.886e+01,
     &  5.946e+01,1.126e+02,5.685e+01,8.710e+01,1.066e+02,1.352e+02,
     &  1.435e+02,1.906e+02,1.908e+02,1.672e+02,1.839e+02,1.708e+02,
     &  1.688e+02,2.351e+02,2.706e+02,3.222e+02,2.620e+02,3.370e+02/
      data (sun3(i),i=1801,1860)/
     &  2.246e+02,2.711e+02,2.740e+02,3.223e+02,2.464e+02,2.250e+02,
     &  2.712e+02,3.674e+02,1.455e+02,1.025e+02,1.822e+02,2.591e+02,
     &  1.733e+02,7.562e+01,1.675e+02,8.840e+01,1.955e+02,9.706e+01,
     &  1.313e+02,1.521e+02,1.425e+02,1.314e+02,2.077e+02,2.148e+02,
     &  2.221e+02,2.609e+02,3.694e+02,2.561e+02,1.989e+02,2.048e+02,
     &  1.751e+02,2.295e+02,1.492e+02,1.115e+02,1.910e+02,2.912e+02,
     &  2.025e+02,2.581e+02,3.009e+02,2.744e+02,2.616e+02,3.026e+02,
     &  2.379e+02,3.989e+02,2.989e+02,2.512e+02,3.752e+02,2.299e+02,
     &  2.515e+02,2.060e+02,2.521e+02,2.668e+02,2.880e+02,2.163e+02,
     &  2.388e+02,2.681e+02,3.094e+02,2.222e+02,2.844e+02,2.888e+02/
      data (sun3(i),i=1861,1920)/
     &  2.518e+02,1.974e+02,2.464e+02,3.416e+02,3.662e+02,1.890e+02,
     &  3.015e+02,2.317e+02,2.568e+02,2.766e+02,2.468e+02,2.684e+02,
     &  2.399e+02,2.336e+02,3.133e+02,2.578e+02,3.153e+02,2.888e+02,
     &  3.633e+02,2.550e+02,2.250e+02,2.327e+02,1.786e+02,2.246e+02,
     &  2.416e+02,2.701e+02,3.772e+02,2.879e+02,2.425e+02,2.990e+02,
     &  2.995e+02,2.105e+02,1.778e+02,1.457e+02,7.720e+01,6.322e+01,
     &  1.064e+02,7.888e+01,1.098e+02,9.072e+01,1.461e+02,1.476e+02,
     &  1.062e+02,1.460e+02,8.945e+01,9.599e+01,1.301e+02,7.194e+01,
     &  7.451e+01,5.351e+01,1.079e+02,1.027e+02,6.458e+01,6.250e+01,
     &  1.209e+02,1.227e+02,1.337e+02,8.561e+01,3.841e+01,5.427e+01/
      data (sun3(i),i=1921,1980)/
     &  1.080e+02,1.280e+02,7.147e+01,1.374e+02,9.922e+01,1.996e+02,
     &  1.352e+02,1.245e+02,5.456e+01,1.015e+02,1.376e+02,1.946e+02,
     &  1.687e+02,1.519e+02,1.291e+02,7.997e+01,8.464e+01,1.539e+02,
     &  1.388e+02,1.282e+02,1.479e+02,1.222e+02,1.136e+02,1.363e+02,
     &  9.164e+01,6.305e+01,6.575e+01,9.692e+01,1.085e+02,9.025e+01,
     &  7.116e+01,7.687e+01,1.037e+02,6.924e+01,9.677e+01,6.201e+01,
     &  4.701e+01,6.277e+01,4.705e+01,6.230e+01,6.935e+01,5.770e+01,
     &  5.294e+01,6.882e+01,5.187e+01,5.422e+01,4.341e+01,5.599e+01,
     &  5.054e+01,6.250e+01,6.156e+01,3.641e+01,3.311e+01,4.122e+01,
     &  4.713e+01,3.622e+01,3.476e+01,4.600e+01,5.285e+01,3.377e+01/
      data (sun3(i),i=1981,2040)/
     &  4.219e+01,3.859e+01,4.106e+01,3.686e+01,5.796e+01,5.901e+01,
     &  5.489e+01,7.432e+01,6.193e+01,3.488e+01,5.671e+01,8.196e+01,
     &  6.374e+01,5.236e+01,5.935e+01,7.155e+01,7.274e+01,6.400e+01,
     &  7.918e+01,6.098e+01,5.309e+01,3.760e+01,2.594e+01,3.493e+01,
     &  3.272e+01,5.058e+01,5.836e+01,3.784e+01,2.624e+01,5.231e+01,
     &  6.284e+01,4.004e+01,5.164e+01,5.762e+01,7.071e+01,6.756e+01,
     &  4.732e+01,5.337e+01,5.979e+01,5.802e+01,5.745e+01,5.218e+01,
     &  5.144e+01,4.906e+01,4.635e+01,3.534e+01,4.767e+01,5.670e+01,
     &  5.252e+01,5.597e+01,4.499e+01,5.471e+01,4.781e+01,4.579e+01,
     &  4.603e+01,5.262e+01,5.848e+01,4.638e+01,6.333e+01,6.969e+01/
      data (sun3(i),i=2041,2100)/
     &  6.198e+01,6.684e+01,8.475e+01,7.516e+01,5.147e+01,6.639e+01,
     &  7.242e+01,5.904e+01,4.354e+01,7.080e+01,5.808e+01,7.827e+01,
     &  7.732e+01,6.127e+01,7.943e+01,8.996e+01,6.309e+01,5.738e+01,
     &  6.335e+01,8.789e+01,7.135e+01,6.003e+01,6.743e+01,6.353e+01,
     &  5.726e+01,5.364e+01,3.072e+01,4.746e+01,2.748e+01,3.586e+01,
     &  4.882e+01,4.172e+01,2.564e+01,2.906e+01,2.866e+01,4.954e+01,
     &  4.829e+01,5.449e+01,3.629e+01,4.069e+01,4.651e+01,2.370e+01,
     &  2.414e+01,4.737e+01,5.656e+01,4.891e+01,6.514e+01,3.573e+01,
     &  3.176e+01,5.154e+01,5.273e+01,4.244e+01,2.847e+01,2.383e+01,
     &  2.527e+01,4.103e+01,3.772e+01,5.802e+01,6.462e+01,4.367e+01/
      data (sun3(i),i=2101,2160)/
     &  4.600e+01,2.743e+01,6.466e+01,6.665e+01,5.024e+01,4.406e+01,
     &  4.300e+01,3.468e+01,4.870e+01,3.738e+01,5.671e+01,4.959e+01,
     &  4.490e+01,2.841e+01,3.720e+01,6.337e+01,5.811e+01,5.665e+01,
     &  4.692e+01,5.812e+01,6.171e+01,5.185e+01,5.224e+01,3.369e+01,
     &  2.443e+01,3.862e+01,2.890e+01,2.843e+01,2.034e+01,4.568e+01,
     &  4.942e+01,4.139e+01,3.106e+01,2.490e+01,4.091e+01,6.734e+01,
     &  4.449e+01,4.418e+01,3.171e+01,3.274e+01,4.878e+01,6.330e+01,
     &  5.296e+01,3.797e+01,2.738e+01,5.176e+01,8.211e+01,7.092e+01,
     &  3.145e+01,3.651e+01,3.577e+01,5.826e+01,4.686e+01,5.272e+01,
     &  6.938e+01,4.323e+01,3.563e+01,3.304e+01,3.057e+01,4.194e+01/
      data (sun3(i),i=2161,2220)/
     &  5.544e+01,5.358e+01,7.675e+01,7.967e+01,5.509e+01,3.162e+01,
     &  5.621e+01,4.262e+01,3.411e+01,3.832e+01,2.636e+01,2.769e+01,
     &  6.848e+01,4.867e+01,4.911e+01,4.705e+01,4.805e+01,4.111e+01,
     &  3.913e+01,5.750e+01,3.711e+01,3.647e+01,7.327e+01,5.108e+01,
     &  4.743e+01,6.879e+01,5.846e+01,4.600e+01,4.756e+01,4.743e+01,
     &  4.324e+01,3.099e+01,4.147e+01,3.100e+01,3.869e+01,3.353e+01,
     &  3.266e+01,3.007e+01,2.612e+01,3.186e+01,3.489e+01,3.122e+01,
     &  2.619e+01,3.952e+01,4.272e+01,5.374e+01,3.969e+01,3.272e+01,
     &  6.732e+01,4.805e+01,5.042e+01,3.893e+01,6.223e+01,3.210e+01,
     &  6.969e+01,4.771e+01,4.838e+01,5.313e+01,5.035e+01,8.284e+01/
      data (sun3(i),i=2221,2280)/
     &  5.962e+01,7.108e+01,4.481e+01,4.740e+01,5.227e+01,4.850e+01,
     &  5.727e+01,7.848e+01,4.580e+01,6.122e+01,6.765e+01,8.359e+01,
     &  6.884e+01,8.798e+01,7.355e+01,4.449e+01,4.375e+01,5.330e+01,
     &  5.738e+01,3.991e+01,3.993e+01,5.581e+01,4.410e+01,5.360e+01,
     &  4.338e+01,5.720e+01,5.419e+01,5.738e+01,3.930e+01,2.626e+01,
     &  1.729e+01,3.990e+01,4.838e+01,4.419e+01,5.060e+01,2.202e+01,
     &  1.686e+01,3.075e+01,4.232e+01,3.349e+01,4.126e+01,5.598e+01,
     &  5.419e+01,6.771e+01,6.643e+01,5.781e+01,2.666e+01,3.649e+01,
     &  5.739e+01,6.306e+01,5.101e+01,3.100e+01,3.587e+01,6.027e+01,
     &  6.500e+01,4.575e+01,3.671e+01,4.481e+01,3.958e+01,4.664e+01/
      data (sun3(i),i=2281,2340)/
     &  3.204e+01,2.037e+01,6.077e+01,4.514e+01,3.080e+01,6.122e+01,
     &  5.890e+01,5.913e+01,5.740e+01,4.575e+01,2.203e+01,2.816e+01,
     &  3.919e+01,2.405e+01,2.529e+01,3.516e+01,4.197e+01,3.179e+01,
     &  4.408e+01,2.210e+01,3.246e+01,2.951e+01,1.504e+01,1.979e+01,
     &  3.906e+01,3.498e+01,4.719e+01,3.388e+01,4.383e+01,4.237e+01,
     &  2.460e+01,2.003e+01,2.041e+01,3.133e+01,3.428e+01,3.913e+01,
     &  5.984e+01,3.579e+01,3.764e+01,1.886e+01,5.297e+01,4.162e+01,
     &  3.732e+01,3.322e+01,3.610e+01,3.482e+01,4.446e+01,5.240e+01,
     &  4.404e+01,4.902e+01,5.670e+01,2.677e+01,1.921e+01,3.275e+01,
     &  4.919e+01,3.033e+01,2.480e+01,2.635e+01,3.751e+01,2.508e+01/
      data (sun3(i),i=2341,2400)/
     &  3.481e+01,2.401e+01,2.568e+01,3.745e+01,3.819e+01,3.457e+01,
     &  3.925e+01,2.334e+01,1.154e+01,1.498e+01,2.897e+01,2.974e+01,
     &  4.410e+01,5.718e+01,3.462e+01,3.537e+01,4.820e+01,4.615e+01,
     &  2.739e+01,2.411e+01,3.083e+01,3.577e+01,4.155e+01,2.980e+01,
     &  1.823e+01,1.893e+01,2.293e+01,3.506e+01,2.732e+01,4.131e+01,
     &  3.547e+01,1.718e+01,1.806e+01,3.905e+01,2.913e+01,1.955e+01,
     &  2.049e+01,2.506e+01,2.768e+01,1.744e+01,2.600e+01,1.984e+01,
     &  1.572e+01,2.821e+01,2.124e+01,1.929e+01,1.685e+01,1.679e+01,
     &  1.421e+01,2.265e+01,1.402e+01,1.903e+01,1.463e+01,1.278e+01,
     &  1.428e+01,1.396e+01,1.367e+01,1.307e+01,1.378e+01,1.407e+01/
      data (sun3(i),i=2401,2460)/
     &  1.296e+01,1.356e+01,1.480e+01,1.382e+01,1.323e+01,1.381e+01,
     &  1.313e+01,1.231e+01,1.190e+01,1.277e+01,1.166e+01,1.112e+01,
     &  1.217e+01,1.218e+01,1.171e+01,1.057e+01,1.253e+01,1.129e+01,
     &  1.238e+01,8.687e+00,1.062e+01,1.212e+01,1.270e+01,1.028e+01,
     &  1.004e+01,1.203e+01,1.108e+01,9.293e+00,1.157e+01,1.150e+01,
     &  1.198e+01,1.150e+01,1.105e+01,1.174e+01,1.049e+01,1.140e+01,
     &  9.622e+00,1.142e+01,1.230e+01,1.116e+01,1.043e+01,1.001e+01,
     &  1.170e+01,1.139e+01,1.078e+01,9.622e+00,8.734e+00,9.449e+00,
     &  1.113e+01,9.840e+00,1.103e+01,9.832e+00,8.365e+00,7.702e+00,
     &  9.123e+00,1.067e+01,1.097e+01,9.437e+00,9.154e+00,8.843e+00/
      data (sun3(i),i=2461,2494)/
     &  1.021e+01,8.275e+00,6.018e+00,7.045e+00,8.009e+00,1.007e+01,
     &  8.980e+00,7.486e+00,7.783e+00,9.095e+00,9.285e+00,8.224e+00,
     &  1.030e+01,8.923e+00,8.154e+00,8.251e+00,8.283e+00,8.517e+00,
     &  9.428e+00,8.940e+00,7.486e+00,7.695e+00,9.173e+00,8.089e+00,
     &  8.193e+00,7.251e+00,8.160e+00,8.710e+00,8.561e+00,7.154e+00,
     &  8.829e+00,8.004e+00,8.322e+00,8.312e+00/

c sun3    modtran3      100-49960 cm-1      20 cm-1           2494

      nns=2494
      wnmin=100.
      wnmax=49960.
      do i=1,nns
        sun(i)=sun3(nns-i+1)
        wn=wnmin+(wnmax-wnmin)*real(nns-i)/(nns-1)
        wlsun(i)=10000./wn
      enddo
      return
      end


c============================================================

      subroutine snow(wlalb,alb,nna)
      use params, only: kr
      implicit none
      integer, parameter :: mxwl=751
      integer :: i, nna
      real(kr) :: albx(mxwl), wlalb(*), alb(*), wlmax, wlmin

c  
c  Wiscombe, W.J. and and S.G. Warren, 1980: "A model for the spectral
c  albedo of snow. I: pure snow."  J. Atmospheric Sciences, 37, 2712-2733.
c  (based on there Mie/delta-eddington model with Reff=125um )

      data (albx(i),i=1,110) /
     & .970,.970,.970,.970,.970,.970,.970,.970,.970,.970,.970,
     & .971,.972,.973,.974,.975,.975,.975,.975,.976,.977,.978,
     & .979,.980,.980,.980,.980,.980,.980,.980,.980,.981,.982,
     & .983,.984,.985,.985,.985,.985,.986,.987,.988,.989,.990,
     & .990,.990,.990,.990,.990,.990,.990,.990,.990,.990,.990,
     & .988,.987,.986,.985,.983,.982,.981,.980,.978,.977,.976,
     & .975,.973,.972,.971,.970,.968,.967,.966,.965,.965,.965,
     & .965,.965,.963,.962,.961,.960,.953,.952,.951,.950,.950,
     & .950,.950,.950,.948,.945,.943,.940,.938,.935,.933,.930,
     & .928,.925,.923,.920,.918,.915,.913,.910,.908,.907,.906/
      data (albx(i),i=111,220) /
     & .905,.903,.902,.901,.900,.898,.895,.893,.890,.885,.880,
     & .875,.870,.865,.860,.855,.850,.843,.842,.841,.840,.835,
     & .830,.825,.820,.810,.800,.790,.780,.771,.762,.753,.745,
     & .739,.733,.727,.720,.718,.715,.713,.710,.715,.720,.725,
     & .730,.733,.735,.738,.740,.742,.743,.744,.745,.745,.745,
     & .745,.745,.743,.742,.741,.740,.732,.724,.716,.710,.692,
     & .674,.656,.640,.630,.620,.610,.600,.584,.568,.552,.535,
     & .524,.513,.502,.490,.490,.490,.490,.490,.490,.490,.490,
     & .490,.492,.493,.494,.495,.495,.495,.495,.495,.495,.495,
     & .495,.495,.493,.492,.491,.490,.488,.485,.483,.480,.470/
      data (albx(i),i=221,330) /
     & .460,.450,.440,.424,.408,.392,.375,.320,.264,.198,.150,
     & .134,.118,.101,.085,.074,.063,.051,.040,.036,.032,.028,
     & .025,.025,.025,.025,.025,.028,.030,.033,.035,.039,.043,
     & .047,.050,.053,.055,.058,.060,.063,.065,.068,.070,.073,
     & .075,.078,.080,.083,.085,.088,.090,.093,.095,.098,.100,
     & .104,.108,.112,.115,.119,.123,.126,.130,.135,.140,.145,
     & .150,.153,.155,.158,.160,.163,.165,.168,.170,.173,.175,
     & .178,.180,.184,.188,.192,.195,.204,.213,.222,.230,.238,
     & .246,.254,.260,.254,.248,.242,.235,.206,.177,.148,.120,
     & .104,.088,.072,.055,.049,.043,.037,.030,.038,.035,.033/
      data (albx(i),i=331,mxwl) /
     & .030,.028,.025,.023,.020,.023,.022,.021,.020,.018,.017,
     & .016,.015,.013,.012,.011,.010,.010,.010,.010,.010,.010,
     & .010,.010,.010,.012,.013,.014,.015,.017,.018,.019,.020,
     & .022,.023,.024,.025,.025,.025,.025,.025,.027,.030,.032,
     & .035,.040,.045,.050,.055,.059,.063,.067,.070,.074,.078,
     & .082,.085,.089,.093,.097,.100,.102,.103,.104,.105,.109,
     & .113,.117,.120,.124,.128,.132,.135,.130,.125,.120,.115,
     & .108,.100,.092,.085,.081,.077,.073,.070,.068,.065,.063,
     & .062,.058,.055,.053,.050,.048,.047,.046,.045,.042,.040,
     & .038,.035,.033,.032,.031,.030,.030,.030,.030,.030,.028,
     & .027,.026,.025,.025,.025,.025,.025,.025,.025,.025,.025,
     & .025,.025,.025,.025,.025,.025,.025,.025,.025,.025,.025,
     & .025,.025,.025,.025,.025,.023,.022,.021,.020,.018,.017,
     & .016,.015,.013,.012,.011,.010,.008,.007,.006,.005,.003,
     & .002,.001,.000,264*0.0/

      nna=mxwl
      alb(1:nna)=albx(1:nna)
      wlmin=.25
      wlmax=4.0
      do i=1,nna
        wlalb(i)=wlmin+(wlmax-wlmin)*real(i-1)/(nna-1)      
      enddo
      return
      end

c=============================================================================

      subroutine clearw(wlalb,alb,nna)
      use params, only: kr
      implicit none
      integer, parameter :: mxwl=751
      integer :: i, nna
      real(kr) :: albx(mxwl), wlalb(*), alb(*), wlmax, wlmin

c     clear water reflectance
c     warning : values of dry sand ground reflectance are given
c     between 0.5 and 1.0 microns. outside this interval the
c     values are set to 0.
c
c     Viollier, M. "Teledetection des concentrations de seston et
c     pigments chlorophylliens contenus dans l'Ocean", These de Doctorat
c     d'Etat, no 503, 1980

      data (albx(i),i=1,mxwl) / 30*0.,
     & .041,.041,.041,.041,.041,.041,.041,.041,.041,.041,.041,
     & .041,.041,.041,.041,.041,.041,.041,.041,.041,.041,.042,
     & .043,.044,.044,.046,.047,.049,.050,.052,.054,.055,.056,
     & .059,.060,.061,.061,.059,.057,.054,.053,.051,.050,.049,
     & .047,.046,.046,.045,.044,.043,.043,.041,.040,.038,.037,
     & .037,.036,.033,.032,.031,.029,.027,.024,.023,.021,.018,
     & .015,.012,.009,.008,.006,.004,.002,.001,.000,.000,.000,
     & 644*0./
      nna=mxwl
      alb(1:nna)=albx(1:nna)
      wlmin=.25
      wlmax=4.0
      do i=1,nna
        wlalb(i)=wlmin+(wlmax-wlmin)*real(i-1)/(nna-1)      
      enddo
      return
      end

c=============================================================================

      subroutine lakew(wlalb,alb,nna)
      use params, only: kr
      implicit none
      integer, parameter :: mxwl=751
      integer :: i, nna
      real(kr) :: albx(mxwl), wlalb(*), alb(*), wlmax, wlmin

c     lake water reflectance
c     warning : values of lake water reflectance are given
c     between 0.35 and 1.0 microns. outside this interval the
c     values are set to 0.
c
c     Kondratyev, K. Y., 1969: "Radiation in the atmosphere", Academic
c     Press, N.Y. 10003 USA


      data (albx(i),i=1,mxwl) / 21*0.,
     & .045,.047,.048,.050,.051,.051,.053,.053,.055,.057,.058,
     & .059,.060,.060,.062,.064,.065,.067,.068,.070,.070,.069,
     & .071,.071,.072,.074,.074,.075,.076,.076,.077,.077,.077,
     & .077,.079,.080,.080,.081,.081,.082,.082,.082,.083,.082,
     & .083,.082,.082,.082,.082,.082,.080,.079,.079,.077,.075,
     & .074,.072,.071,.070,.068,.066,.065,.064,.063,.061,.060,
     & .059,.057,.056,.054,.053,.051,.050,.048,.047,.046,.045,
     & .044,.042,.041,.040,.039,.037,.036,.035,.034,.032,.032,
     & .031,.030,.030,.028,.027,.027,.025,.024,.024,.023,.022,
     & .022,.020,.020,.020,.019,.019,.019,.019,.019,.019,.018,
     & .016,.016,.015,.014,.012,.012,.012,.011,.011,.009,.008,
     & .008,.006,.005,.005,.003,.002,.000,602*0./
      nna=mxwl
      alb(1:nna)=albx(1:nna)
      wlmin=.25
      wlmax=4.0
      do i=1,nna
        wlalb(i)=wlmin+(wlmax-wlmin)*real(i-1)/(nna-1)      
      enddo
      return
      end

c=============================================================================

      subroutine seaw(wlalb,alb,nna)
      use params, only: kr
      implicit none
      integer, parameter :: mxwl=751
      integer :: i, nna
      real(kr) :: albx(mxwl), wlalb(*), alb(*), wlmax, wlmin
c
c  sea water (actually just clear water)
c
c     Viollier, M. "Teledetection des concentrations de seston et
c     pigments chlorophylliens contenus dans l'Ocean", These de Doctorat
c     d'Etat, no 503, 1980


      data (albx(i),i=1,mxwl) / 41*0.041,
     & .041,.041,.041,.041,.041,.041,.041,.041,.041,.041,.042,
     & .043,.044,.044,.046,.047,.049,.050,.052,.054,.055,.056,
     & .059,.060,.061,.061,.059,.057,.054,.053,.051,.050,.049,
     & .047,.046,.046,.045,.044,.043,.043,.041,.040,.038,.037,
     & .037,.036,.033,.032,.031,.029,.027,.024,.023,.021,.018,
     & .015,.012,.009,.008,.006,.004,.002,.001,.000,.000,.000,
     & 644*0./
      nna=mxwl
      alb(1:nna)=albx(1:nna)
      wlmin=.25
      wlmax=4.0
      do i=1,nna
        wlalb(i)=wlmin+(wlmax-wlmin)*real(i-1)/(nna-1)      
      enddo
      return
      end

c=============================================================================

      subroutine sand(wlalb,alb,nna)
      use params, only: kr
      implicit none
      integer, parameter :: mxwl=751
      integer :: i, nna
      real(kr) :: albx(mxwl), wlalb(*), alb(*), wlmax, wlmin

c     sand average reflectance
c     warning : values of dry sand ground reflectance are given
c     between 0.4 and 2.2 microns. outside this interval the
c     values are set to 0.
c
c     Staetter, R., and M. Schroeder, 1978: "Spectral characteristics of
c     natural surfaces", Proceeding of the tenth Int. Conf. on Earth
c     Obs. from Space, 6-11 March 1978, (ESA-SP, 134)


      data (albx(i),i=1,140) / 30*0.,
     & .091,.091,.091,.091,.091,.091,.095,.095,.095,.095,.095,
     & .095,.095,.097,.097,.100,.103,.103,.107,.107,.107,.107,
     & .110,.114,.117,.121,.121,.125,.125,.125,.128,.126,.131,
     & .134,.134,.134,.134,.134,.134,.141,.145,.149,.149,.154,
     & .157,.160,.163,.163,.166,.169,.173,.177,.181,.185,.189,
     & .195,.199,.204,.204,.208,.214,.220,.224,.228,.233,.235,
     & .239,.242,.246,.246,.248,.251,.254,.256,.260,.263,.263,
     & .263,.266,.270,.273,.276,.279,.279,.282,.284,.284,.284,
     & .286,.286,.286,.290,.292,.292,.292,.292,.292,.292,.292,
     & .295,.295,.298,.298,.301,.301,.306,.306,.310,.310,.314/
      data (albx(i),i=141,240) /
     & .314,.314,.316,.316,.316,.316,.316,.319,.321,.321,
     & .321,.321,.324,.324,.324,.324,.322,.322,.322,.320,
     & .320,.317,.317,.317,.320,.320,.320,.323,.323,.327,
     & .327,.327,.329,.329,.332,.332,.332,.332,.333,.335,
     & .338,.341,.345,.348,.353,.353,.353,.353,.355,.355,
     & .355,.355,.353,.353,.353,.353,.353,.353,.356,.356,
     & .356,.356,.360,.360,.360,.363,.363,.363,.367,.367,
     & .367,.370,.370,.370,.370,.370,.370,.367,.367,.364,
     & .362,.362,.362,.359,.359,.359,.359,.359,.363,.363,
     & .363,.363,.363,.365,.369,.369,.369,.369,.369,.372/
      data (albx(i),i=241,340) /
     & .372,.372,.372,.375,.375,.375,.375,.375,.375,.375,
     & .375,.375,.375,.375,.375,.375,.375,.375,.377,.377,
     & .377,.377,.377,.377,.377,.377,.377,.377,.377,.377,
     & .379,.379,.379,.379,.382,.382,.382,.382,.385,.385,
     & .385,.385,.385,.385,.388,.388,.388,.388,.388,.388,
     & .388,.391,.391,.391,.391,.391,.391,.391,.391,.391,
     & .391,.391,.394,.394,.394,.394,.394,.394,.394,.396,
     & .396,.396,.396,.396,.396,.393,.393,.393,.393,.389,
     & .386,.386,.382,.382,.378,.378,.374,.369,.369,.369,
     & .371,.371,.371,.371,.371,.371,.371,.371,.371,.371/
      data (albx(i),i=341,mxwl) /
     & .374,.374,.377,.377,.377,.377,.379,.379,.379,.379,
     & .383,.385,.388,.388,.388,.392,.392,.392,.395,.395,
     & .395,.395,.395,.393,.393,.393,.388,.388,.388,.388,
     & .388,.385,.385,.385,.381,.381,.381,.381,.381,.381,
     & .381,.381,.374,.374,.374,.374,.374,.374,.374,.374,
     & .372,.369,.369,.369,.369,.369,.369,.369,.369,.369,
     & .369,.369,.369,.369,.369,.369,.369,.369,.369,.369,
     & .369,340*0./
      nna=mxwl
      alb(1:nna)=albx(1:nna)
      wlmin=.25
      wlmax=4.0
      do i=1,nna
        wlalb(i)=wlmin+(wlmax-wlmin)*real(i-1)/(nna-1)      
      enddo
      return
      end

c=============================================================================

      subroutine vegeta(wlalb,alb,nna)
      use params, only: kr
      implicit none
      integer, parameter :: mxwl=751
      integer :: i, nna
      real(kr) :: albx(mxwl), wlalb(*), alb(*), wlmax, wlmin

c     vegetation average reflectance
c     warning : values of dry sand ground reflectance are given
c     between 0.4 and 2.2 microns. outside this interval the
c     values are set to 0.
c
c     Manual of Remote Sensing American Society of Photogrammetry,
c     R.G. Reeves, A. Anson, D. Landen, eds.  1st ed. Falls Church, Va.,
c     1975.

      data (albx(i),i=1,140) / 30*0.,
     & .060,.060,.062,.064,.066,.069,.070,.072,.074,.078,.080,
     & .083,.084,.089,.093,.098,.102,.104,.106,.110,.115,.119,
     & .121,.125,.127,.130,.133,.133,.134,.133,.131,.127,.121,
     & .115,.110,.105,.101,.098,.094,.090,.087,.083,.081,.080,
     & .078,.076,.075,.074,.073,.073,.073,.074,.079,.100,.138,
     & .169,.199,.228,.259,.290,.316,.350,.378,.403,.436,.462,
     & .487,.509,.511,.514,.519,.520,.520,.522,.522,.522,.523,
     & .524,.524,.524,.524,.526,.526,.526,.527,.527,.527,.528,
     & .528,.528,.529,.529,.529,.529,.531,.531,.531,.531,.531,
     & .532,.532,.532,.532,.532,.533,.533,.533,.534,.534,.534/
      data (albx(i),i=141,240) /
     & .534,.535,.535,.536,.536,.537,.537,.536,.536,.535,
     & .535,.534,.532,.531,.530,.528,.528,.527,.527,.526,
     & .525,.524,.522,.521,.519,.518,.515,.513,.512,.510,
     & .508,.507,.506,.505,.502,.500,.498,.496,.495,.493,
     & .492,.492,.492,.492,.492,.493,.495,.495,.496,.496,
     & .496,.497,.497,.497,.498,.498,.497,.497,.497,.495,
     & .493,.492,.491,.488,.486,.482,.478,.476,.472,.467,
     & .462,.451,.441,.429,.421,.408,.399,.385,.371,.365,
     & .349,.339,.330,.321,.309,.298,.289,.279,.272,.267,
     & .259,.251,.243,.233,.229,.224,.218,.215,.215,.215/
      data (albx(i),i=241,340) /
     & .215,.219,.223,.229,.234,.240,.249,.256,.260,.267,
     & .273,.279,.286,.293,.300,.306,.312,.319,.325,.331,
     & .337,.341,.345,.351,.355,.360,.362,.367,.369,.372,
     & .376,.378,.379,.381,.382,.384,.386,.387,.389,.388,
     & .388,.388,.388,.388,.388,.384,.383,.381,.380,.378,
     & .376,.374,.373,.371,.370,.368,.367,.366,.365,.365,
     & .363,.362,.361,.359,.358,.357,.355,.353,.350,.347,
     & .346,.345,.343,.340,.337,.335,.331,.330,.321,.312,
     & .296,.273,.221,.186,.158,.138,.129,.121,.110,.102,
     & .095,.091,.089,.086,.086,.084,.084,.084,.086,.087/
      data (albx(i),i=341,mxwl) /
     & .093,.098,.105,.114,.116,.124,.133,.134,.141,.147,
     & .151,.156,.162,.166,.170,.174,.175,.178,.181,.185,
     & .187,.188,.192,.196,.199,.201,.205,.208,.212,.213,
     & .214,.217,.219,.220,.221,.224,.227,.229,.231,.233,
     & .237,.238,.239,.241,.242,.243,.245,.245,.246,.248,
     & .248,.250,.246,.242,.238,.234,.230,.226,.222,.218,
     & .214,.210,.206,.202,.198,.194,.190,.186,.182,.178,
     & .174,.170,.166,.162,.158,.154,.150,.146,.142,.138,
     & .134,.130,.126,.122,.118,.114,.110,.106,.102,.098,
     & .094,.090,.086,.082,.078,.074,.070,.066,.062,.058,
     & .054,.050,.046,.042,.038,.034,.030,.026,.022,.018,
     & .014,.010,.010,.010,.010,.010,.010,.010,.010,.010,
     & .010,.010,.010,.010,.010,.010,.010,.010,.010,.010,
     & .010,280*0./
      nna=mxwl
      alb(1:nna)=albx(1:nna)
      wlmin=.25
      wlmax=4.0
      do i=1,nna
        wlalb(i)=wlmin+(wlmax-wlmin)*real(i-1)/(nna-1)      
      enddo
      return
      end

c=======================================================================
c=======================================================================
      module fltblk             ! filter data base
      use params, only: kr
      save
      integer, parameter :: mxwl=1000
      real(kr), dimension(mxwl) :: wlfilt,filt
      integer :: nnf=0
      end module fltblk
c=======================================================================
      subroutine setfilt(isat,wlinf,wlsup,wlinc,wlmin,wlmax,nwl)
c
c input   isat      filter selection flag
c 
c         wlinf     lower wavelength limit
c                      or central wavelength for isat=-2,-3
c         wlsup     upper wavelength limit
c                      or equivalent width for isat=-2,-3
c         wlinc     wavelength increment, used to determine number
c                      of wavelength sample points in wavelength iteration
c output
c         wlmin     lower wavelength limit
c         wlmax     upper wavelength limit
c         nwl       number of sample points
c
c         wlinc     is set to 0.001 if single wavelength sample point 
c                     is specified, i.e., either if
c
c                     isat.eq.-2 .and. wlsup.eq.0   
c                          -or-
c                     isal.eq.0 .and. wlinf.eq.wlsup
c
c                   otherwise it is not changed
c
c NOTE: actual filter values are returned by the function "filter"
c
      use params, only: kr, pi
      use fltblk, only: nnf, mxwl, filt, wlfilt
      implicit none
      integer :: i, isat, nwl
      real(kr) :: sqpi, wlinc, wlinf, wlmax, wlmin, wlsup, xc, xlim, xx 

      if(wlsup.eq.0.and.isat.lt.-2) then
        write(*,*) 'Error -- WLSUP must be non-zero when ISAT=',isat
        write(*,*) '         WLSUP =',wlsup
        stop
      endif

      select case (isat)

      case (-4)                 ! center, equivalent width, gaussian
        sqpi=sqrt(pi)
        xc=2
        xlim=xc*sqpi
        nnf=1000
        wlmin=wlinf-xc*wlsup
        wlmax=wlinf+xc*wlsup
        do i=1,nnf
          xx=-xlim+(i-1)*(2*xlim)/(nnf-1)
          filt(i)=exp(-xx**2)
          wlfilt(i)=wlmin+(wlmax-wlmin)*real(i-1)/(nnf-1)
        enddo

      case (-3)               ! center, equivalent width, triangular
        nnf=3
        wlmin=wlinf-wlsup
        wlmax=wlinf+wlsup
        filt(1:3)=(/0.,1.,0./)
        wlfilt(1:3)=(/wlmin,wlinf,wlmax/)

      case(-2)               ! center, width, flat top
        nnf=0
        wlmin=wlinf-.5*wlsup
        wlmax=wlinf+.5*wlsup
        if(wlsup.eq.0.) then
          nwl=1
          wlinc=.001
          return
        endif

      case (-1)              ! read from filter.dat
        nnf=mxwl
        call rdspec('filter.dat',wlfilt,filt,nnf)
        wlmin=wlfilt(1)
        wlmax=wlfilt(nnf)

      case (0)               ! wl_min, wl_max, flat  top
        nnf=0
        wlmin=wlinf
        wlmax=wlsup
        if(wlinf.eq.wlsup) then
          nwl=1
          wlinc=.001
          return
        endif

      case (1)  ;  call meteo(filt,wlmin,wlmax,nnf)
      case (2)  ;  call goese(filt,wlmin,wlmax,nnf)
      case (3)  ;  call goesw(filt,wlmin,wlmax,nnf)
      case (4)  ;  call avhr81(filt,wlmin,wlmax,nnf)
      case (5)  ;  call avhr82(filt,wlmin,wlmax,nnf)
      case (6)  ;  call avhr91(filt,wlmin,wlmax,nnf)
      case (7)  ;  call avhr92(filt,wlmin,wlmax,nnf)
      case (8)  ;  call avhr101(filt,wlmin,wlmax,nnf)
      case (9)  ;  call avhr102(filt,wlmin,wlmax,nnf)
      case (10) ;  call avhr111(filt,wlmin,wlmax,nnf)
      case (11) ;  call avhr112(filt,wlmin,wlmax,nnf)
      case (12) ;  call gtr1(filt,wlmin,wlmax,nnf)
      case (13) ;  call gtr2(filt,wlmin,wlmax,nnf)
      case (14) ;  call nm410(filt,wlmin,wlmax,nnf)
      case (15) ;  call nm936(filt,wlmin,wlmax,nnf)
      case (16) ;  call mfrsr1(filt,wlmin,wlmax,nnf)
      case (17) ;  call mfrsr2(filt,wlmin,wlmax,nnf)
      case (18) ;  call mfrsr3(filt,wlmin,wlmax,nnf)
      case (19) ;  call mfrsr4(filt,wlmin,wlmax,nnf)
      case (20) ;  call mfrsr5(filt,wlmin,wlmax,nnf)
      case (21) ;  call mfrsr6(filt,wlmin,wlmax,nnf)
      case (22) ;  call avhr83(filt,wlmin,wlmax,nnf)
      case (23) ;  call avhr84(filt,wlmin,wlmax,nnf)
      case (24) ;  call avhr85(filt,wlmin,wlmax,nnf)
      case (25) ;  call setlow(filt,wlmin,wlmax,nnf)
      case (26) ;  call airs1(filt,wlmin,wlmax,nnf)
      case (27) ;  call airs2(filt,wlmin,wlmax,nnf)
      case (28) ;  call airs3(filt,wlmin,wlmax,nnf)
      case (29) ;  call airs4(filt,wlmin,wlmax,nnf)

      end select

      if(wlmin.lt.0.199) then
        write(*,*) 'Error in SETFILT -- illegal wavelength limits '
        write(*,*) wlmin,wlmax
        stop
      endif

      if(isat.gt.0) then
        do i=1,nnf
          wlfilt(i)=wlmin+(wlmax-wlmin)*real(i-1)/(nnf-1)
        enddo
      endif
      
      if(wlinc.gt.1.) then
c       equal spacing in wavenumber
        nwl=((10000./wlmin) - (10000./wlmax))/wlinc + 1.
      elseif(wlinc.lt.0.) then
c       equal spacing in ln(wavelength) and ln(wavenumber)
        nwl=1+log(wlmax/wlmin)/abs(wlinc)
      else
c       use default wavelength increment
        if(wlinc.eq.0)
     &       wlinc=(wlmax-wlmin)/max(10,1+int((wlmax-wlmin)/0.005))
c       equal spacing in wavelength
        nwl=nint((wlmax-wlmin)/wlinc)+1
      endif

      if(wlmin.ne.wlmax.and.nwl.eq.1) nwl=2

      return
      end     
      
c=======================================================================
      function filter(w)
c
c input:  w       wavelength in microns
c output:         filter functions
c 
c NOTE: "setfilt" must be called before "filter"
c
      use params, only: one, zero, kr
      use fltblk, only: nnf, filt, wlfilt
      implicit none
      integer :: i
      real(kr) :: w, wt, filter

      if(nnf.eq.0) then 
        filter=1.
      else 
        call locate(wlfilt,nnf,w,i)
        wt=(w-wlfilt(i))/(wlfilt(i+1)-wlfilt(i))
        wt=max(zero,min(one,wt))
        filter=filt(i)*(1.-wt)+filt(i+1)*wt
      endif
      return
      end
c=======================================================================
      subroutine meteo(srr,wmin,wmax,nnf)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=.355, wmx=1.11
      integer, parameter :: mwv=152
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax
c
c    spectral response of meteosat
c
      data sr/ 
     & 0.0050,0.0094,0.0154,0.0191,0.0208,0.0263,0.0357,0.0404,0.0513,
     & 0.0570,0.0689,0.0746,0.0820,0.0900,0.1001,0.1071,0.1152,0.1222,
     & 0.1358,0.1482,0.1561,0.1683,0.1789,0.1908,0.2047,0.2146,0.2280,
     & 0.2426,0.2582,0.2758,0.2897,0.3029,0.3209,0.3388,0.3514,0.3705,
     & 0.3886,0.4035,0.4240,0.4387,0.4644,0.4885,0.5095,0.5348,0.5611,
     & 0.5841,0.6109,0.6332,0.6478,0.6617,0.6736,0.6857,0.7006,0.7120,
     & 0.7254,0.7430,0.7698,0.7834,0.8002,0.8206,0.8434,0.8563,0.8753,
     & 0.8949,0.9214,0.9410,0.9571,0.9722,0.9844,0.9893,0.9948,0.9968,
     & 0.9960,1.0000,0.9960,0.9953,0.9921,0.9874,0.9824,0.9770,0.9735,
     & 0.9710,0.9660,0.9660,0.9623,0.9623,0.9628,0.9599,0.9569,0.9509,
     & 0.9405,0.9271,0.9150,0.9029,0.8865,0.8753,0.8644,0.8548,0.8426,
     & 0.8283,0.8151,0.8037,0.7901,0.7804,0.7660,0.7591,0.7489,0.7393,
     & 0.7326,0.7192,0.6984,0.6664,0.6434,0.6178,0.6005,0.5824,0.5603,
     & 0.5385,0.5232,0.5038,0.4862,0.4706,0.4471,0.4211,0.4007,0.3809,
     & 0.3618,0.3388,0.3147,0.2924,0.2758,0.2548,0.2332,0.2099,0.1908,
     & 0.1678,0.1509,0.1348,0.1138,0.0994,0.0838,0.0751,0.0642,0.0538,
     & 0.0456,0.0372,0.0263,0.0178,0.0097,0.0084,0.0020,0.0000/
c
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
      return
      end
c=======================================================================
      subroutine goese(srr,wmin,wmax,nnf)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=.500, wmx=.895
      integer, parameter :: mwv=80
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax
c
c    spectral response of goes east
c
      data sr/
     & 0.013,0.046,0.133,0.302,0.455,0.584,0.680,0.727,0.772,0.807,
     & 0.838,0.859,0.879,0.899,0.919,0.934,0.948,0.965,0.983,0.991,
     & 0.983,0.974,0.961,0.949,0.937,0.925,0.916,0.912,0.907,0.911,
     & 0.914,0.907,0.892,0.877,0.843,0.808,0.772,0.734,0.695,0.678,
     & 0.660,0.639,0.613,0.588,0.556,0.521,0.486,0.457,0.427,0.397,
     & 0.363,0.330,0.301,0.279,0.257,0.234,0.210,0.186,0.165,0.146,
     & 0.127,0.110,0.096,0.082,0.069,0.058,0.047,0.038,0.031,0.024,
     & 0.018,0.014,0.010,0.006,0.004,0.003,0.002,0.001,0.001,0.000/
c
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
      return
      end
c=======================================================================
      subroutine goesw(srr,wmin,wmax,nnf)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=.500, wmx=.905
      integer, parameter :: mwv=82
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax
c
c    spectral response of goes west
c
      data sr/
     & 0.002,0.020,0.065,0.183,0.409,0.611,0.775,0.889,0.925,0.957,
     & 0.973,0.987,0.993,0.999,0.998,0.997,0.990,0.983,0.974,0.964,
     & 0.956,0.951,0.944,0.928,0.912,0.893,0.874,0.853,0.832,0.811,
     & 0.795,0.780,0.760,0.738,0.716,0.686,0.656,0.625,0.592,0.560,
     & 0.536,0.513,0.489,0.464,0.440,0.413,0.385,0.358,0.331,0.305,
     & 0.279,0.256,0.233,0.211,0.192,0.174,0.157,0.143,0.129,0.115,
     & 0.101,0.087,0.077,0.070,0.062,0.054,0.046,0.038,0.032,0.027,
     & 0.022,0.018,0.015,0.011,0.008,0.006,0.004,0.003,0.002,0.001,
     & 0.001,0.000/
c
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
      return
      end
c=======================================================================
      subroutine avhr81(srr,wmin,wmax,nnf)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=0.505, wmx=0.745
      integer, parameter :: mwv=49
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax
c
c    1st spectral response of avhrr  (noaa 8)
c
      data sr/
     &  .004,.008,.019,.028,.035,.048, .060, .070, .078, .087,
     &  .095,.099,.367,.527,.665,.788, .815, .811, .807, .828,
     &  .870,.903,.931,.956,.978,.993,1.000,1.000,1.000,1.000,
     & 1.000,.976,.900,.817,.743,.622, .438, .306, .217, .140,
     &  .088,.067,.047,.032,.024,.015, .003, .001, .000/
c
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c=======================================================================
      subroutine avhr82(srr,wmin,wmax,nnf)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=0.70, wmx=1.085
      integer, parameter :: mwv=78
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax
c
c    2nd spectral response of avhr  (noaa 8)
c
      data sr/
     & .035,.168,.350,.510,.622,.756,.837,.896,.929,.949,
     & .952,.946,.945,.986,.998,.998,.998,.988,.939,.910,
     & .889,.869,.853,.840,.827,.817,.809,.802,.797,.793,
     & .789,.785,.781,.775,.768,.759,.753,.749,.745,.741,
     & .738,.719,.695,.678,.662,.648,.636,.627,.618,.612,
     & .603,.600,.597,.593,.587,.582,.562,.523,.485,.433,
     & .376,.288,.230,.189,.145,.107,.075,.051,.030,.020,
     & .014,.008,.006,.003,.001,.001,.001,.000/
c
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c=======================================================================
      subroutine avhr91(srr,wmin,wmax,nnf)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=.54, wmx=.80
      integer, parameter :: mwv=53
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax
c
c     revised  jan94
c
c    1st spectral response of avhr (noaa 9)
c
      data sr/
     &    0.0010, 0.0040, 0.0400, 0.1398, 0.2700, 0.3929, 0.5000,
     &    0.5898, 0.6700, 0.7448, 0.8000, 0.8203, 0.8150, 0.7973,
     &    0.7850, 0.7933, 0.8200, 0.8587, 0.8950, 0.9134, 0.9050,
     &    0.8663, 0.8250, 0.8087, 0.8150, 0.8397, 0.8900, 0.9645,
     &    1.0000, 0.9363, 0.8100, 0.6633, 0.5000, 0.3232, 0.1800,
     &    0.1138, 0.0900, 0.0697, 0.0500, 0.0351, 0.0250, 0.0188,
     &    0.0150, 0.0124, 0.0100, 0.0073, 0.0050, 0.0038, 0.0030,
     &    0.0019, 0.0010, 0.0009, 0.0010/

c
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c=======================================================================
      subroutine avhr92(srr,wmin,wmax,nnf)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=.645, wmx=1.190
      integer, parameter :: mwv=110
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax
c
c    2nd spectral response of avhr  (noaa 9)
c
      data sr/
     & .0036,.0071,.0071,.0071,.0075,.0079,.0079,.0079,.0100,.0120,
     & .0864,.1607,.3715,.4222,.5381,.6540,.7888,.9235,.9481,.9727,
     & .9845,.9963,.9982,1.000,.9945,.9891,.9891,.9891,.9593,.9294,
     & .9366,.9438,.9276,.9114,.8919,.8724,.8724,.8724,.8724,.8724,
     & .8724,.8724,.8724,.8724,.8724,.8724,.8710,.8695,.8647,.8600,
     & .8533,.8465,.8418,.8372,.8284,.8197,.8041,.7886,.7845,.7805,
     & .7762,.7719,.7656,.7594,.7374,.7154,.6489,.5824,.4963,.4103,
     & .3259,.2415,.1862,.1309,.1009,.0709,.0537,.0364,.0299,.0234,
     & .0188,.0142,.0106,.0070,.0070,.0070,.0061,.0051,.0038,.0025,
     & .0025,.0025,.0025,.0025,.0025,.0025,.0025,.0025,.0025,.0025,
     & .0025,.0025,.0025,.0025,.0025,.0025,.0025,.0025,.0025,.0000/
c
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c=======================================================================
      subroutine avhr101(srr,wmin,wmax,nnf)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=.555, wmx=.750
      integer, parameter :: mwv=40
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data sr/ 
     $    0.039998,    0.079998,    0.274987,    0.469986,    0.604990,
     $    0.739990,    0.779997,    0.819997,    0.815000,    0.810000,
     $    0.829999,    0.849999,    0.894997,    0.939997,    0.964998,
     $    0.989998,    0.975001,    0.960001,    0.955000,    0.950000,
     $    0.974998,    0.999998,    0.995000,    0.990000,    0.810015,
     $    0.630015,    0.455016,    0.280017,    0.192508,    0.105008,
     $    0.077503,    0.050003,    0.040001,    0.030001,    0.025000,
     $    0.020000,    0.015000,    0.010000,    0.008000,    0.000000/
c
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
       return
       end
c======================================================================       
      subroutine avhr102(srr,wmin,wmax,nnf)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=.685, wmx=1.180
      integer, parameter :: mwv=100
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data sr(1:50)/
     $    0.001000,    0.002000,    0.002000,    0.002000,    0.040996,
     $    0.079996,    0.239985,    0.399985,    0.564984,    0.729984,
     $    0.804993,    0.879993,    0.909997,    0.939997,    0.959998,
     $    0.979998,    0.989999,    0.999999,    0.985002,    0.970002,
     $    0.955002,    0.940002,    0.920002,    0.900002,    0.894501,
     $    0.889001,    0.883501,    0.878001,    0.872501,    0.867001,
     $    0.861501,    0.856001,    0.850501,    0.845001,    0.838001,
     $    0.831001,    0.824001,    0.817001,    0.810001,    0.803001,
     $    0.796001,    0.789001,    0.782001,    0.775001,    0.763501,
     $    0.752001,    0.740501,    0.729001,    0.717502,    0.706002/
      data sr(51:100)/
     $    0.694502,    0.683002,    0.671502,    0.660002,    0.617006,
     $    0.574006,    0.531006,    0.488006,    0.445006,    0.402006,
     $    0.359006,    0.316007,    0.273006,    0.230006,    0.211503,
     $    0.193003,    0.174503,    0.156003,    0.137503,    0.119003,
     $    0.100003,    0.081003,    0.062503,    0.044003,    0.025503,
     $    0.007003,    0.006500,    0.006000,    0.005500,    0.005000,
     $    0.004000,    0.003000,    0.002500,    0.002000,    0.002000,
     $    0.002000,    0.002000,    0.002000,    0.001500,    0.001000,
     $    0.001000,    0.001000,    0.001000,    0.001000,    0.001000,
     $    0.001000,    0.001000,    0.001000,    0.000500,    0.000000/
c
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c======================================================================       
      subroutine avhr111(srr,wmin,wmax,nnf)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=.545, wmx=.750
      integer, parameter :: mwv=42
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data sr/
     $    0.000500,    0.001000,    0.047997,    0.094997,    0.297487,
     $    0.499986,    0.569995,    0.639995,    0.699996,    0.759996,
     $    0.779999,    0.799999,    0.785001,    0.770001,    0.779999,
     $    0.789999,    0.839996,    0.889996,    0.885000,    0.880000,
     $    0.850002,    0.820002,    0.802501,    0.785002,    0.799999,
     $    0.814999,    0.907492,    0.999992,    0.900009,    0.800010,
     $    0.605019,    0.410019,    0.290011,    0.170011,    0.122505,
     $    0.075005,    0.057502,    0.040002,    0.030001,    0.020001,
     $    0.017000,    0.000000/
c
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
       return
       end
c======================================================================       
      subroutine avhr112(srr,wmin,wmax,nnf)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=.685, wmx=1.100
      integer, parameter :: mwv=84
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data sr/
     $    0.004000,    0.008000,    0.018999,    0.029999,    0.189985,
     $    0.349985,    0.494986,    0.639986,    0.749990,    0.859990,
     $    0.899996,    0.939996,    0.964997,    0.989997,    0.994000,
     $    0.998000,    0.983002,    0.968002,    0.946502,    0.925002,
     $    0.892503,    0.860004,    0.855001,    0.850001,    0.847500,
     $    0.845000,    0.842500,    0.840000,    0.837500,    0.835000,
     $    0.832500,    0.830000,    0.827500,    0.825000,    0.824500,
     $    0.824000,    0.823500,    0.823000,    0.822500,    0.822000,
     $    0.821500,    0.821000,    0.820500,    0.820000,    0.801003,
     $    0.782003,    0.763003,    0.744003,    0.725003,    0.706003,
     $    0.687003,    0.668003,    0.649003,    0.630003,    0.585006,
     $    0.540006,    0.495006,    0.450006,    0.405006,    0.360006,
     $    0.315007,    0.270007,    0.225007,    0.180006,    0.165502,
     $    0.151002,    0.136002,    0.121003,    0.106502,    0.092002,
     $    0.077502,    0.063002,    0.048003,    0.033003,    0.018502,
     $    0.004002,    0.003500,    0.003000,    0.002500,    0.002000,
     $    0.001500,    0.001000,    0.000500,    0.000000/
c
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c======================================================================       
      subroutine gtr1(srr,wmin,wmax,nnf)
c      gtr-100 ch 1 (bsi's predeployment calibration measurement 
c      linearly interpolated to 5nm intervals)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=.555, wmx=.805
      integer, parameter :: mwv=51
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data sr/
     &    0.0048,    0.0501,    0.2077,    0.4466,    0.6530,
     &    0.7900,    0.8780,    0.9232,    0.9468,    0.9712,
     &    0.9877,    0.9930,    0.9944,    1.0000,    0.9940,
     &    0.9785,    0.9700,    0.9538,    0.9241,    0.9531,
     &    0.9024,    0.8780,    0.8550,    0.8311,    0.8164,
     &    0.7743,    0.7516,    0.7424,    0.6844,    0.6400,
     &    0.5675,    0.5034,    0.4515,    0.3408,    0.2231,
     &    0.1521,    0.1057,    0.0716,    0.0511,    0.0414,
     &    0.0332,    0.0247,    0.0320,    0.0226,    0.0195,
     &    0.0167,    0.0153,    0.0139,    0.0113,    0.0117,
     &    0.0000/
c
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c======================================================================       
      subroutine gtr2(srr,wmin,wmax,nnf)
c      gtr-100 ch 2 (sbrc's postdeployment calibration measurements
c      spline fit to 5nm intervals)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=.705, wmx=1.130
      integer, parameter :: mwv=86
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data sr(1:50)/
     &   0.0713, 0.1737, 0.2703, 0.3611, 0.4461,
     &   0.5248, 0.5972, 0.6629, 0.7213, 0.7720,
     &   0.8146, 0.8500, 0.8790, 0.9026, 0.9217,
     &   0.9368, 0.9490, 0.9586, 0.9666, 0.9734,
     &   0.9797, 0.9854, 0.9904, 0.9946, 0.9977,
     &   0.9996, 1.0000, 0.9988, 0.9956, 0.9903,
     &   0.9828, 0.9735, 0.9628, 0.9511, 0.9387,
     &   0.9261, 0.9134, 0.9012, 0.8897, 0.8793,
     &   0.8704, 0.8626, 0.8555, 0.8489, 0.8424,
     &   0.8356, 0.8281, 0.8196, 0.8096, 0.7978/
      data sr(51:86)/
     &   0.7839, 0.7679, 0.7499, 0.7300, 0.7084,
     &   0.6851, 0.6603, 0.6338, 0.6057, 0.5760,
     &   0.5447, 0.5122, 0.4787, 0.4447, 0.4106,
     &   0.3766, 0.3432, 0.3106, 0.2792, 0.2493,
     &   0.2212, 0.1950, 0.1704, 0.1476, 0.1263,
     &   0.1067, 0.0886, 0.0721, 0.0572, 0.0440,
     &   0.0325, 0.0227, 0.0144, 0.0077, 0.0025,
     &   0.0000/
c
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c======================================================================       
      subroutine nm410(srr,wmin,wmax,nnf)
c      gtr-100 410nm channel
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=.389, wmx=.423   
      integer, parameter :: mwv=35
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data sr/
     &   0.0000,    0.0124,    0.0116,    0.0157,    0.0101,    0.0293,
     &   0.0280,    0.0510,    0.0549,    0.0566,    0.0797,    0.1107,
     &   0.1277,    0.1723,    0.2512,    0.3532,    0.5287,    0.7189,
     &   0.8874,    0.9670,    1.0000,    0.9932,    0.9853,    0.9357,
     &   0.8522,    0.6905,    0.4765,    0.2993,    0.1549,    0.0775,
     &   0.0400,    0.0263,    0.0131,    0.0035,    0.0000/
c
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c=======================================================================
      subroutine nm936(srr,wmin,wmax,nnf)
c      gtr-100 410nm channel
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=.878, wmx=.955
      integer, parameter :: mwv=78
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data sr/
     &  0.0004,    0.0018,    0.0038,    0.0059,    0.0071,    0.0070,
     &  0.0070,    0.0080,    0.0095,    0.0117,    0.0134,    0.0149,
     &  0.0152,    0.0160,    0.0187,    0.0246,    0.0324,    0.0389,
     &  0.0383,    0.0349,    0.0383,    0.0416,    0.0383,    0.0373,
     &  0.0436,    0.0621,    0.0888,    0.1028,    0.0833,    0.0734,
     &  0.0849,    0.0926,    0.0928,    0.0824,    0.0775,    0.0888,
     &  0.1136,    0.1431,    0.1580,    0.1482,    0.1292,    0.1255,
     &  0.1252,    0.1264,    0.1263,    0.1236,    0.1229,    0.1270,
     &  0.1392,    0.1690,    0.2288,    0.3384,    0.4835,    0.6716,
     &  0.8228,    0.9022,    0.9113,    0.9200,    0.9597,    1.0000,
     &  0.9954,    0.9095,    0.7721,    0.5937,    0.4003,    0.2264,
     &  0.1207,    0.0613,    0.0323,    0.0172,    0.0100,    0.0062,
     &  0.0042,    0.0030,    0.0022,    0.0014,    0.0010,    0.0000/
c
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c=======================================================================
      subroutine mfrsr1(srr,wmin,wmax,nnf)
c      MFRSR channel 1
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=0.397, wmx=0.433
      integer, parameter :: mwv=37
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data sr/
     &    0.0001,    0.0003,    0.0008,    0.0020,    0.0044,
     &    0.0092,    0.0185,    0.0349,    0.0625,    0.1058,
     &    0.1696,    0.2570,    0.3686,    0.5000,    0.6417,
     &    0.7792,    0.8950,    0.9727,    1.0000,    0.9727,
     &    0.8950,    0.7792,    0.6417,    0.5000,    0.3686,
     &    0.2570,    0.1696,    0.1058,    0.0625,    0.0349,
     &    0.0185,    0.0092,    0.0044,    0.0020,    0.0008,
     &    0.0003,    0.0001/
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c=======================================================================
      subroutine mfrsr2(srr,wmin,wmax,nnf)
c      MFRSR channel 2
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=0.482, wmx=0.518
      integer, parameter :: mwv=37
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data sr/
     &    0.0001,    0.0003,    0.0008,    0.0020,    0.0044,
     &    0.0092,    0.0185,    0.0349,    0.0625,    0.1058,
     &    0.1696,    0.2570,    0.3686,    0.5000,    0.6417,
     &    0.7792,    0.8950,    0.9727,    1.0000,    0.9727,
     &    0.8950,    0.7792,    0.6417,    0.5000,    0.3686,
     &    0.2570,    0.1696,    0.1058,    0.0625,    0.0349,
     &    0.0185,    0.0092,    0.0044,    0.0020,    0.0008,
     &    0.0003,    0.0001/
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c=======================================================================
      subroutine mfrsr3(srr,wmin,wmax,nnf)
c      MFRSR channel 3
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=0.592, wmx=0.628
      integer, parameter :: mwv=37
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data sr/
     &    0.0001,    0.0003,    0.0008,    0.0020,    0.0044,
     &    0.0092,    0.0185,    0.0349,    0.0625,    0.1058,
     &    0.1696,    0.2570,    0.3686,    0.5000,    0.6417,
     &    0.7792,    0.8950,    0.9727,    1.0000,    0.9727,
     &    0.8950,    0.7792,    0.6417,    0.5000,    0.3686,
     &    0.2570,    0.1696,    0.1058,    0.0625,    0.0349,
     &    0.0185,    0.0092,    0.0044,    0.0020,    0.0008,
     &    0.0003,    0.0001/
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c=======================================================================
      subroutine mfrsr4(srr,wmin,wmax,nnf)
c      MFRSR channel 4
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=0.647, wmx=0.683
      integer, parameter :: mwv=37
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data sr/
     &    0.0001,    0.0003,    0.0008,    0.0020,    0.0044,
     &    0.0092,    0.0185,    0.0349,    0.0625,    0.1058,
     &    0.1696,    0.2570,    0.3686,    0.5000,    0.6417,
     &    0.7792,    0.8950,    0.9727,    1.0000,    0.9727,
     &    0.8950,    0.7792,    0.6417,    0.5000,    0.3686,
     &    0.2570,    0.1696,    0.1058,    0.0625,    0.0349,
     &    0.0185,    0.0092,    0.0044,    0.0020,    0.0008,
     &    0.0003,    0.0001/
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c=======================================================================
      subroutine mfrsr5(srr,wmin,wmax,nnf)
c      MFRSR channel 5
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=0.844, wmx=0.880
      integer, parameter :: mwv=37
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data sr/
     &    0.0001,    0.0003,    0.0008,    0.0020,    0.0044,
     &    0.0092,    0.0185,    0.0349,    0.0625,    0.1058,
     &    0.1696,    0.2570,    0.3686,    0.5000,    0.6417,
     &    0.7792,    0.8950,    0.9727,    1.0000,    0.9727,
     &    0.8950,    0.7792,    0.6417,    0.5000,    0.3686,
     &    0.2570,    0.1696,    0.1058,    0.0625,    0.0349,
     &    0.0185,    0.0092,    0.0044,    0.0020,    0.0008,
     &    0.0003,    0.0001/
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c=======================================================================
      subroutine mfrsr6(srr,wmin,wmax,nnf)
c      MFRSR channel 6
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=0.922, wmx=0.958
      integer, parameter :: mwv=37
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data sr/
     &    0.0001,    0.0003,    0.0008,    0.0020,    0.0044,
     &    0.0092,    0.0185,    0.0349,    0.0625,    0.1058,
     &    0.1696,    0.2570,    0.3686,    0.5000,    0.6417,
     &    0.7792,    0.8950,    0.9727,    1.0000,    0.9727,
     &    0.8950,    0.7792,    0.6417,    0.5000,    0.3686,
     &    0.2570,    0.1696,    0.1058,    0.0625,    0.0349,
     &    0.0185,    0.0092,    0.0044,    0.0020,    0.0008,
     &    0.0003,    0.0001/
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c=======================================================================
      subroutine avhr83(srr,wmin,wmax,nnf)
c     noaa-8    ahvrr channel 3
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=3.40, wmx=4.05
      integer, parameter :: mwv=27
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data sr/
     &   0.0000, 0.0084, 0.0953, 0.2715, 0.4908, 0.7088, 0.8807, 0.9713,
     &   1.0000, 0.9932, 0.9719, 0.9551, 0.9506, 0.9427, 0.9295, 0.9197,
     &   0.9128, 0.9067, 0.8991, 0.8869, 0.8325, 0.6741, 0.3732, 0.1003,
     &   0.0306, 0.0056, 0.0000/
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c=======================================================================
      subroutine avhr84(srr,wmin,wmax,nnf)
c     noaa-8    ahvrr channel 4
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=10.225, wmx=11.675
      integer, parameter :: mwv=30
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data sr/
     &   0.0000, 0.1228, 0.3305, 0.5884, 0.8287, 0.9459, 0.9889, 1.0000,
     &   0.9901, 0.9670, 0.9371, 0.9061, 0.8765, 0.8479, 0.8193, 0.7899,
     &   0.7590, 0.7255, 0.6881, 0.6450, 0.5895, 0.5109, 0.3960, 0.2395,
     &   0.1127, 0.0509, 0.0238, 0.0101, 0.0013, 0.0000/
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c=======================================================================
      subroutine avhr85(srr,wmin,wmax,nnf)
c     noaa-8    ahvrr channel 5
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=11.300, wmx=12.600
      integer, parameter :: mwv=27
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data sr/
     &   0.0000, 0.0613, 0.2334, 0.5287, 0.8555, 0.9537, 0.9249, 0.9351,
     &   0.9946, 1.0000, 0.9491, 0.8988, 0.8862, 0.8969, 0.9098, 0.9071,
     &   0.8867, 0.8515, 0.8029, 0.7409, 0.6642, 0.5720, 0.4648, 0.3412,
     &   0.2026, 0.0694, 0.0000/
      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
c
      return
      end
c=======================================================================

      subroutine setlow(srr,wmin,wmax,nnf)
c     action spectra for DNA damage by UVB radiation (Setlow, R.B, 1974)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=.280, wmx=0.341
      integer, parameter :: mwv=61

      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data  sr /
     $  7.6459E-01, 6.8615E-01, 6.1575E-01, 5.5258E-01, 4.9589E-01,
     $  4.4501E-01, 3.9935E-01, 3.5838E-01, 3.2161E-01, 2.8862E-01,
     $  2.5901E-01, 2.1868E-01, 1.8459E-01, 1.5582E-01, 1.3153E-01,
     $  1.1103E-01, 8.7078E-02, 6.8320E-02, 5.3602E-02, 4.2055E-02,
     $  3.2996E-02, 2.3459E-02, 1.6681E-02, 1.1862E-02, 8.4351E-03,
     $  5.9981E-03, 4.1937E-03, 2.9307E-03, 2.0480E-03, 1.4312E-03,
     $  1.0002E-03, 6.9896E-04, 4.8845E-04, 3.4134E-04, 2.3854E-04,
     $  1.6670E-04, 1.1649E-04, 8.1409E-05, 5.6891E-05, 3.9757E-05,
     $  2.7783E-05, 1.9416E-05, 1.3568E-05, 9.4818E-06, 6.6261E-06,
     $  4.6305E-06, 3.2359E-06, 2.2614E-06, 1.5803E-06, 1.1044E-06,
     $  7.7176E-07, 5.3932E-07, 3.7689E-07, 2.6338E-07, 1.8406E-07,
     $  1.2863E-07, 8.9888E-08, 6.2816E-08, 4.3898E-08, 3.0677E-08,
     $  2.1438E-08 /

      wmin=wmn
      wmax=wmx
      nnf=mwv
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
      return
      end
c=======================================================================
C spectral response of AIRS channel 1
C=======================================================================
      subroutine airs1(srr,wmin,wmax,nnf)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=0.380, wmx=0.460
      integer,  parameter :: mwv= 35
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data (sr  (i),i=  1, 35) /
     $  0.005128, 0.002886, 0.001374, 0.003698, 0.006864,
     $  0.013090, 0.033231, 0.076092, 0.161704, 0.245375,
     $  0.259199, 0.267823, 0.274717, 0.286149, 0.298826,
     $  0.309475, 0.317824, 0.317501, 0.310042, 0.307051,
     $  0.292431, 0.294837, 0.338277, 0.363464, 0.327798,
     $  0.239252, 0.149011, 0.085511, 0.046337, 0.024886,
     $  0.014469, 0.009347, 0.007380, 0.005831, 0.004215 /

       wmin=wmn
       wmax=wmx
       nnf=mwv

       do 10 i=1,nnf
         srr(i)=sr(i)
 10   continue
       return
       end

C=======================================================================
C spectral response of AIRS channel 2
C=======================================================================
      subroutine airs2(srr,wmin,wmax,nnf)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=0.520, wmx=0.700
      integer,  parameter :: mwv= 75
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax

      data (sr  (i),i=    1,   50) / 
     $  0.000000, 0.000533, 0.000945, 0.001361, 0.001783,
     $  0.002211, 0.002644, 0.003083, 0.003527, 0.003976,
     $  0.004142, 0.004254, 0.004366, 0.004958, 0.007938,
     $  0.010861, 0.018275, 0.027527, 0.047194, 0.072460,
     $  0.112121, 0.173171, 0.264447, 0.374261, 0.510745,
     $  0.652248, 0.738061, 0.773845, 0.777518, 0.778838,
     $  0.784899, 0.792303, 0.796726, 0.793400, 0.789225,
     $  0.789243, 0.791918, 0.794664, 0.793183, 0.782235,
     $  0.774498, 0.763391, 0.764142, 0.767907, 0.782163,
     $  0.794344, 0.804426, 0.810698, 0.811115, 0.809902 /

      data (sr  (i),i= 51, 75) / 
     $  0.815325, 0.826189, 0.836451, 0.852932, 0.852806,
     $  0.844644, 0.833519, 0.828699, 0.841141, 0.851687,
     $  0.851055, 0.846426, 0.825965, 0.713892, 0.432026,
     $  0.216361, 0.114826, 0.063494, 0.038689, 0.021043,
     $  0.011658, 0.006644, 0.004261, 0.001918, 0.000000 /
       wmin=wmn
       wmax=wmx
       nnf=mwv

       do 10 i=1,nnf
         srr(i)=sr(i)
 10   continue
       return
       end

C=======================================================================
C spectral response of AIRS channel 3
C=======================================================================
      subroutine airs3(srr,wmin,wmax,nnf)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=0.670, wmx=0.975
      integer,  parameter :: mwv=125
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax
      
      data (sr  (i),i=    1,   50) / 
     $     0.000000, 0.003379, 0.003411, 0.003444, 0.003477,
     $     0.004078, 0.006951, 0.011480, 0.016599, 0.023960,
     $     0.036043, 0.057332, 0.088698, 0.148032, 0.220794,
     $     0.325841, 0.444009, 0.547685, 0.660726, 0.775603,
     $     0.853849, 0.906058, 0.918434, 0.916672, 0.907973,
     $     0.901806, 0.908660, 0.918596, 0.928232, 0.935273,
     $     0.937799, 0.940109, 0.942381, 0.944376, 0.947208,
     $     0.956354, 0.963515, 0.969355, 0.973644, 0.977545,
     $     0.974620, 0.968519, 0.960392, 0.952298, 0.945408,
     $     0.938639, 0.932385, 0.930490, 0.928068, 0.923987 /
      
      data (sr  (i),i=   51,  100) / 
     $     0.917406, 0.908386, 0.898606, 0.890075, 0.880939,
     $     0.875101, 0.873721, 0.875880, 0.879767, 0.885273,
     $     0.890469, 0.894408, 0.887393, 0.879201, 0.866884,
     $     0.850730, 0.835656, 0.824285, 0.815378, 0.807416,
     $     0.804449, 0.803807, 0.804861, 0.805892, 0.807177,
     $     0.807003, 0.804723, 0.795822, 0.784649, 0.776168,
     $     0.765594, 0.752811, 0.744292, 0.737139, 0.735197,
     $     0.731744, 0.726148, 0.715571, 0.707805, 0.702938,
     $     0.697320, 0.691714, 0.684362, 0.676938, 0.667892,
     $     0.656787, 0.648466, 0.641009, 0.634255, 0.626658 /
      
      data (sr  (i),i=101,125) / 
     $     0.620257, 0.611595, 0.600326, 0.591389, 0.580594,
     $     0.565131, 0.546977, 0.530590, 0.514660, 0.504639,
     $     0.502949, 0.509015, 0.509009, 0.469097, 0.392224,
     $     0.264232, 0.148339, 0.070440, 0.035891, 0.016855,
     $     0.009114, 0.005118, 0.003445, 0.002607, 0.000000 /
      wmin=wmn
      wmax=wmx
      nnf=mwv
      
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
      return
      end
      
C=======================================================================
C spectral response of AIRS channel 4
C=======================================================================
      subroutine airs4(srr,wmin,wmax,nnf)
      use params, only: kr
      implicit none
      real(kr), parameter :: wmn=0.415, wmx=1.110
      integer,  parameter :: mwv=259
      integer :: i,nnf
      real(kr) :: srr(*),sr(mwv),wmin,wmax
      
      data (sr  (i),i=    1,   50) / 
     $     0.000000, 0.000835, 0.001808, 0.002801, 0.003818,
     $     0.004859, 0.008237, 0.014668, 0.028408, 0.053174,
     $     0.226002, 0.461725, 0.491061, 0.458212, 0.422702,
     $     0.417594, 0.427578, 0.449909, 0.483016, 0.512304,
     $     0.533860, 0.546248, 0.556111, 0.554623, 0.553001,
     $     0.554188, 0.555123, 0.562048, 0.570187, 0.581164,
     $     0.593501, 0.603637, 0.611861, 0.620140, 0.628474,
     $     0.636258, 0.639973, 0.643660, 0.647318, 0.650875,
     $     0.657180, 0.664420, 0.671698, 0.679299, 0.688091,
     $     0.696974, 0.705916, 0.713593, 0.720423, 0.727282 /
      
      data (sr  (i),i=   51,  100) / 
     $     0.734167, 0.736732, 0.739088, 0.741418, 0.743194,
     $     0.743984, 0.744725, 0.745639, 0.746794, 0.747931,
     $     0.749048, 0.750125, 0.753838, 0.762325, 0.770839,
     $     0.778745, 0.785975, 0.792448, 0.798838, 0.805255,
     $     0.811696, 0.818163, 0.819922, 0.821466, 0.822827,
     $     0.821900, 0.820401, 0.818876, 0.816191, 0.812171,
     $     0.808835, 0.805440, 0.805978, 0.807495, 0.808469,
     $     0.809383, 0.810284, 0.811171, 0.813025, 0.817150,
     $     0.821234, 0.825328, 0.829432, 0.837509, 0.845904,
     $     0.854338, 0.862465, 0.870015, 0.876968, 0.883365 /
      
      data (sr  (i),i=  101,  150) / 
     $     0.889779, 0.896213, 0.901125, 0.904839, 0.908558,
     $     0.912286, 0.913276, 0.913860, 0.914441, 0.914464,
     $     0.913303, 0.912134, 0.911242, 0.910351, 0.909405,
     $     0.908314, 0.906574, 0.904577, 0.903073, 0.901843,
     $     0.900614, 0.899385, 0.898151, 0.896913, 0.895868,
     $     0.897468, 0.899067, 0.900667, 0.902265, 0.904212,
     $     0.906268, 0.908322, 0.910378, 0.912433, 0.915588,
     $     0.919173, 0.922757, 0.926341, 0.929925, 0.933505,
     $     0.934687, 0.935659, 0.936632, 0.936368, 0.935805,
     $     0.935232, 0.934783, 0.934669, 0.934546, 0.933097 /
      
      data (sr  (i),i=  151,  200) / 
     $     0.930430, 0.927784, 0.925141, 0.922344, 0.918725,
     $     0.913528, 0.908044, 0.902574, 0.897117, 0.892657,
     $     0.888827, 0.885007, 0.880852, 0.875619, 0.870432,
     $     0.865259, 0.860098, 0.854950, 0.849800, 0.844656,
     $     0.839525, 0.832753, 0.826716, 0.822859, 0.818949,
     $     0.814976, 0.809735, 0.804009, 0.798282, 0.792555,
     $     0.787259, 0.782719, 0.778161, 0.773586, 0.768983,
     $     0.764229, 0.759209, 0.754179, 0.743430, 0.734616,
     $     0.729730, 0.724829, 0.719912, 0.713143, 0.706277,
     $     0.698287, 0.688177, 0.680948, 0.674610, 0.667855 /
      
      data (sr  (i),i=  201,  250) / 
     $     0.659749, 0.651586, 0.641477, 0.630669, 0.622755,
     $     0.614839, 0.604736, 0.592914, 0.583120, 0.573417,
     $     0.561795, 0.547182, 0.534179, 0.521942, 0.495328,
     $     0.483003, 0.470702, 0.458425, 0.446218, 0.433655,
     $     0.419071, 0.405505, 0.392962, 0.370931, 0.358746,
     $     0.346591, 0.334465, 0.314971, 0.303614, 0.295342,
     $     0.285874, 0.273863, 0.254465, 0.243492, 0.235679,
     $     0.205100, 0.183823, 0.166805, 0.147939, 0.131689,
     $     0.119560, 0.108231, 0.097669, 0.087686, 0.078087,
     $     0.071138, 0.064213, 0.058244, 0.052492, 0.046142 /
      
      data (sr  (i),i=251,259) / 
     $     0.039644, 0.035735, 0.031886, 0.027426, 0.025916,
     $     0.022519, 0.018541, 0.006526, 0.000000 /
      wmin=wmn
      wmax=wmx
      nnf=mwv
      
      do 10 i=1,nnf
        srr(i)=sr(i)
 10   continue
      return
      end
c=======================================================================
      subroutine rdspec(file,wl,r,nn)
c
c purpose:     read user specified spectral input from files
c              filter.dat, albedo.dat and solar.dat
c input:
c  file        input file name
c  nn          maximum number of input values 
c              (an error is triggered if the number of 
c              values read from file exceeds nn)
c
c output:
c  wl          wavelength array
c  r           spectral array
c  nn          number of values read from file
c
      use params, only: kr
      implicit none
      character*(*) file
      character*(80) line
      integer :: nn, nnsv, i,j
      real(kr) :: wl(*),r(*)
c mfl
c     print *, 'in redspec()'     
      wl(1:nn)=0.
      nnsv=nn
      open(unit=13,file=file,status='old',form='formatted')

      read(13,*,end=10) (wl(i),r(i),i=1,huge(0))
 10   continue

      close(unit=13)

      nn=count(wl(1:nn).ne.0.)
      
c mfl 
c     do i=1,nn
c       print *, file, wl(i),r(i)
c     enddo

      if(nn.gt.nnsv) then
c        print *,'ERROR in rdspec --- too many values read from ',file
        print *,'file should not specify more than ',nnsv,' values'
        stop
      endif

      if(wl(1).gt.wl(nn)) then  ! reverse if necessary
c        print *, 'Reversing wavelength array!!!!!!'
        wl(nn:1:-1)=wl(1:nn)
        r(nn:1:-1)=r(1:nn)
      endif

      return
      end
      
c============================================================
      subroutine zensun(iday,time,alat,alon,zenith,azimuth,solfac)
c
c   routine:      zensun
c  
c   purpose:  compute the solar zenith and azimuth angles and solar flux
c             multiplier for a given location, time and day.
c
c   input:
c     iday    day of year (used as a fraction of 365.24 days year)
c  
c     time    universal time in decimal hours
c  
c     alat    geographic latitude of point on earth's surface (degrees)
c  
c     alon    geographic longitude of point on earth's surface (degrees)
c  
c   output:
c  
c     zenith  solar zenith angle (degrees)
c  
c     azimuth solar azimuth measured clockwise from due north (degrees)
c
c     solfac  solar flux multiplier. solfac=1./rsun**2 
c             where rsun is the current earth-sun distance in 
c             astronomical units.  
c  
      use params, only: pi, kr
      implicit none
      
      integer :: nday(74), iday, i

      real(kr), parameter ::
     &     degpday=360./365.242, ! degree per day
     &     eccen=0.01671,        ! eccentricity of earth orbit
     &     dayph=2.              ! day of perihelion


      real(kr) :: eqt(74), dec(74), dtor, frac, eqtime, decang,
     &     sunlat, sunlon, t0, t1, p0, alon, alat, p1, zz, xx, yy,
     &     azimuth, zenith, rsun, solfac, dd, time

      data nday/
     &     1,   6,  11,  16,  21,  26,  31,  36,  41,
     &     46,  51,  56,  61,  66,  71,  76,  81,  86,
     &     91,  96, 101, 106, 111, 116, 121, 126, 131,
     &     136, 141, 146, 151, 156, 161, 166, 171, 176,
     &     181, 186, 191, 196, 201, 206, 211, 216, 221,
     &     226, 231, 236, 241, 246, 251, 256, 261, 266,
     &     271, 276, 281, 286, 291, 296, 301, 306, 311,
     &     316, 321, 326, 331, 336, 341, 346, 351, 356,
     &     361, 366/
      
      data eqt/
     &     -3.23, -5.49, -7.60, -9.48,-11.09,-12.39,-13.34,-13.95,
     &     -14.23,-14.19,-13.85,-13.22,-12.35,-11.26,-10.01, -8.64,
     &     -7.18, -5.67, -4.16, -2.69, -1.29, -0.02,  1.10,  2.05,
     &     2.80,  3.33,  3.63,  3.68,  3.49,  3.09,  2.48,  1.71,
     &     0.79, -0.24, -1.33, -2.41, -3.45, -4.39, -5.20, -5.84,
     &     -6.28, -6.49, -6.44, -6.15, -5.60, -4.82, -3.81, -2.60,
     &     -1.19,  0.36,  2.03,  3.76,  5.54,  7.31,  9.04, 10.69,
     &     12.20, 13.53, 14.65, 15.52, 16.12, 16.41, 16.36, 15.95,
     &     15.19, 14.09, 12.67, 10.93,  8.93,  6.70,  4.32,  1.86,
     &     -0.62, -3.23/
      
      data dec/
     &     -23.06,-22.57,-21.91,-21.06,-20.05,-18.88,-17.57,-16.13,
     &     -14.57,-12.91,-11.16, -9.34, -7.46, -5.54, -3.59, -1.62,
     &     0.36,  2.33,  4.28,  6.19,  8.06,  9.88, 11.62, 13.29,
     &     14.87, 16.34, 17.70, 18.94, 20.04, 21.00, 21.81, 22.47,
     &     22.95, 23.28, 23.43, 23.40, 23.21, 22.85, 22.32, 21.63,
     &     20.79, 19.80, 18.67, 17.42, 16.05, 14.57, 13.00, 11.33,
     &     9.60,  7.80,  5.95,  4.06,  2.13,  0.19, -1.75, -3.69,
     &     -5.62, -7.51, -9.36,-11.16,-12.88,-14.53,-16.07,-17.50,
     &     -18.81,-19.98,-20.99,-21.85,-22.52,-23.02,-23.33,-23.44,
     &     -23.35,-23.06/
      
cccc
c
c compute the subsolar coordinates
c
      dtor=pi/180.
      
      dd=mod(iday-1, 365)+1
      do 10 i=1,74
        if(nday(i) .gt. dd ) goto 20
 10   continue
 20   continue
      frac=(dd-nday(i-1))/(nday(i)-nday(i-1))
      eqtime=eqt(i-1)*(1.-frac)+frac*eqt(i)
      decang=dec(i-1)*(1.-frac)+frac*dec(i)
      sunlat=decang
      sunlon=-15.*(time-12.+eqtime/60.)
      
c      write(12,'(2i4,9es11.3)') i,nday(i),frac,eqtime,sunlat,sunlon
c      write(12,'(9es11.3)') eqt(i-1),eqt(i),eqtime
c      write(12,'(9es11.3)') dec(i-1),dec(i),decang
c
c compute the solar zenith, azimuth and flux multiplier
c
      t0=(90.-alat)*dtor                            
      t1=(90.-sunlat)*dtor                         
      p0=alon*dtor                                  
      p1=sunlon*dtor                               
      zz=cos(t0)*cos(t1)+sin(t0)*sin(t1)*cos(p1-p0) 
c      write(12,'(9es11.3)') t0,t1,p0,p1,zz
c      write(12,'(9es11.3)') cos(t0),cos(t1),sin(t0),sin(t1),cos(p1-p0)
      xx=sin(t1)*sin(p1-p0)
      yy=sin(t0)*cos(t1)-cos(t0)*sin(t1)*cos(p1-p0)
      azimuth=atan2(xx,yy)/dtor
      zenith=acos(zz)/dtor                         
c
      rsun=1.-eccen*cos(degpday*(dd-dayph)*dtor)      
      solfac=1./rsun**2                              
      return
      end
      
      
c=======================================================================

