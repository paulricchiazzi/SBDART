      module outblk
      use params, only: mxly, kr, nstrms

      integer, parameter :: maxulv=mxly+1

      real(kr) :: phidw,topdn,topup,topdir,botdn,botup,botdir,
     &        fxup(maxulv),fxdn(maxulv),fxdir(maxulv),
     &        uurs(nstrms,nstrms),uurl(nstrms,nstrms,maxulv)

      integer, parameter :: nu1=nstrms*nstrms,
     &                      nu2=nu1*maxulv

      data uurs/nu1*0.0/
      data uurl/nu2*0.0/
      data phidw,topdn,topup,topdir,botdn,botup,botdir/7*0./
      data fxup,fxdn,fxdir/maxulv*0.,maxulv*0.,maxulv*0./

      end module outblk

c=======================================================================
c sbdart
c
c PURPOSE: COMPUTE PLANE-PARALLEL RADIATIVE TRANSFER IN THE EARTH'S ATMOSPHERE
c
c INPUT/OUTPUT VARIABLES   described in rt.doc
c
c modules:
c
c   params       grid sizes, physical constants
c   outblk       output quantities
c   albblk       surface albedo quantities
c   sunblk       solar spectra 
c   fltblk       sensor filter function
c   aeroblk      aerosol properties
c   gasblk       gaseous absorption quantities
c   trcblk       gaseous absorption due to trace gases
c   kdstblk      k-distribution quantities
c 
c
c  subroutines called from the main routine
c
c   absint       lowtran gas absorption integrals
c   atms         atmospheric profile
c   bdrefchk     check brdf function
c   chkin        check input parameters
c   chkprn       diagnostic print out
c   depthscl     combine gas absorption with scattering optical depth
c   disort       discrete ordinate radiative tranfer module
c   gasinit      set up k-distribution parameters
c   gasset       set up k-distribution parameters for a particular iteration
c   helper       output format description 
c   locate       search for a value in a 1-d monotonic array
c   modatm       modify atmospheric profile
c   modmix       modify profile of trace gases
c   nearest      find index of array element closest to a given value
c   normom       normalize phase function moments
c   rayleigh     find rayleigh scattering optical depth
c   satcloud     put saturate water vapor in cloudy parts of atmosphere
c   saturate     find fully saturated water vapor density
c   setfilt      initialize sensor filter function
c   stdout0      write head records of standard output format
c   stdout1      write wavelength dependent output records
c   stdout2      write cumulative output record(s)
c   suralb       compute surface albedo 
c   tauaero      compute aerosol optical properties
c   taucloud     compute cloud optical properties
c   usrcloud     read file for cloud properties
c   vuangles     set radiance viewing angles
c   wllimits     coordinate wavelength iteration
c   zensun       solar ephemeris
c   zeroit       zero out an array
c   zgrid        change vertical resolution of atmospheric grid
c   zlayer       find layer positions given atititude ranges
c
c logical units            files
c
c    11           INPUT (stays open)
c
c    12           CKATM and CKTAU k-distribution files (stays open)
c
c    13           usrcld.dat, atms.dat, filter.dat, albedo.dat,
c                 solar.dat (each file is closed after read)
c
c    15           aerosol.dat (stays open)
c
c    16           disort warning messages
c 
c=======================================================================

      program sbdart

      use params, only: mxly, nstrms, ncldz, mxkd, mxq, maxmom,
     &     zip, nan, kr, one, zero, pi

      use outblk, only: maxulv ! passes information to stdout2

      use aeroblk, only:
     &     zaer,taerst,zbaer,dbaer,vis,tbaer,abaer,
     &     wlbaer,qbaer,wbaer,gbaer,pmaer,rhaer,iaer,imoma,jaer

      implicit none

      integer, parameter ::
     &     nta=9,               ! number of gas transmission values
     &     ndb=20               ! idb array size

      integer :: nnn, nz, idatm, isat, nf, isalb, iout, idb(ndb)=0,
     &     imomc, iday, nosct, ntau, numu, nzen, nphi, nothrm, ibcnd,
     &     kdist, ngrid, krhclr, ipth, maxumu, maxphi, nstr, ios, nvzen,
     &     numset, i, ii, j, k, nstrsv, mcldz, nbot, ntop, nwl, il,
     &     nmom, kd, nk, lcld(ncldz), ib=1, nb=1

      real(kr) :: z(mxly), p(mxly), t(mxly), wh(mxly), wo(mxly),
     &     gwk(mxkd), dtauk(mxly,mxkd), dtaugc(mxly), dtauc(mxly),
     &     wcld(mxly), dtaus(mxly), wreal(mxly), dtor, dtaur(mxly),
     &     dtaua(mxly), waer(mxly), dtaug(mxly), umu(nstrms),
     &     utau(mxly), uzen(nstrms), phi(nstrms), vzen(nstrms), sc(5),
     &     uu(mxq, mxly), zout(2),
     &     zcloud(ncldz), tcloud(ncldz), pmom(0:maxmom,mxly),
     &     temper(0:mxly), lwp(ncldz), nre(ncldz), wvnmlo, wvnmhi,
     &     accur, wlinf, wlsup, wlinc, albcon, sza, csza, rhcld, amix,
     &     alat, alon, time, solfac, sclh2o, pbar, uw, uo3, o3trp, ztrp,
     &     xn2, xo2, xco2, xch4, xn2o, xco, xno2, xso2, xnh3, xno,
     &     xhno3, xo4, xrsc, zpres, phi0, fisot, temis, btemp, ttemp,
     &     zgrid1, zgrid2, saza, fj, amu0, wl1, wl2, wl, dw, dwl, ff,
     &     flxin, test, rsfc, tauaertot, wt, dtauctot, rfldir(maxulv),
     &     rfldn(maxulv), flup(maxulv), dfdt(maxulv), uavg(maxulv),
     &     uur(nstrms,maxulv,nstrms), albmed(nstrms), trnmed(nstrms),
     &     ewcoef, etirr, dref, szabar


      logical :: prnt(7), lamber, azmavg, usrtau, usrang, plank,
     &     onlyfl, radcalc, zflag=.false., corint=.false.,
     &     spowder=.false., first=.true.

      real(kr), external :: filter, relhum, solirr, salbedo, albave

      character*(127) header

c=======================================================================

      data wvnmlo/0.0/,wvnmhi/0.0/,accur/0.0/

      data  nz, idatm,isat,nf,isalb,iout,imomc
     &    / 33,   4,   0,  2,   0,   10,   3  /
      data wlinf,   wlsup,   wlinc
     &     / .55,     .55,     0./
      data    albcon,   sza, saza, csza, rhcld, amix
     &     /    0.,    0.,   180.,  zip,   zip,   zip /

      data  dtaua  /mxly*0./
      data  waer   /mxly*0./

      data iday,alat,alon,time      /0,-64.767,-64.067,16./
      data solfac                   /1.0/

      data sclh2o,  pbar,    uw,   uo3, o3trp,  ztrp
     &     /  zip,   zip,   zip,   zip,   zip,   0./
      data    xn2,   xo2,  xco2,  xch4,  xn2o
     &     /  zip,   zip,   zip,   zip,   zip/
      data    xco,  xno2,  xso2,  xnh3,   xno,  xhno3
     &     /  zip,   zip,   zip,   zip,   zip,  zip/

      data xo4/one/
      data xrsc/one/

      data nosct/0/

      data zpres/zip/
c             
      data  zcloud,    tcloud,     nre,     lcld
     &     /ncldz*0., ncldz*0.,  ncldz*8., ncldz*0/

      data prnt/ 7*.false./
      data lamber,  azmavg, usrang, usrtau
     &   / .true.,  .true., .false., .false./

      data ntau,numu,nzen,nphi/4*0/
      data umu/nstrms*0./ 
      data uzen/nstrms*zip/
      data vzen/nstrms*90./
      data phi/nstrms*zip/
      data nothrm,ibcnd,phi0,fisot,temis/-1,0,0.,0.,0./
      data btemp,ttemp /zip,zip/
      data zout/0.,100./ 
      data kdist/3/

      data zgrid1,zgrid2,ngrid          /1.,30.,0/
      data krhclr /0/

      data ipth /0/

      data maxumu,maxphi/2*nstrms/
      data nstr/0/

      data lwp/ncldz*0./

      data sc/5*nan/

      namelist /input/
     &     idatm, amix, isat, wlinf, wlsup, wlinc, sza, csza, solfac,
     &     nf, iday, time, alat, alon, zpres, pbar, sclh2o, uw, uo3,
     &     o3trp, ztrp, xrsc, xn2, xo2, xco2, xch4, xn2o, xco, xno2,
     &     xso2, xnh3, xno, xhno3, xo4, isalb, albcon, sc, zcloud,
     &     tcloud, lwp, nre, rhcld, krhclr, jaer, zaer, taerst, iaer,
     &     vis, rhaer, tbaer, wlbaer, qbaer, abaer, wbaer, gbaer, pmaer,
     &     zbaer, dbaer, nothrm, nosct, kdist, zgrid1, zgrid2, ngrid,
     &     idb, zout, iout, prnt, temis, nstr, nzen, uzen, vzen, nphi,
     &     phi, saza, imomc, imoma, ttemp, btemp, corint, spowder



      namelist /dinput/
     &     ibcnd,phi0,prnt,ipth,fisot,temis,nstr,
     &     nzen,uzen,vzen,nphi,phi,ttemp,btemp


c-----------------------------------------------------------------------
      
      open(unit=11,file='INPUT',status='old',iostat=ios)
      if(ios.eq.0) then
        read(11,input,end=1)
        read(11,dinput,end=2)
        goto 2
 1      continue
        stop 'error: namelist block $INPUT not found'
 2      continue
      else
        write(*,input)
        stop
      endif
      rewind 11
      dtor=pi/180.

      if(idb(1).ne.0) call helper(iout)
 
      radcalc=any(iout.eq.(/5,6,20,21,22,23/))

      onlyfl=.not. radcalc

      if(nstr.eq.0) then
        if(radcalc) then
          nstr=min(20,nstrms)
        else
          nstr=4
        endif
      endif

      select case (abs(isalb))
      case  (7)
        ! chlor, wndspd, salin
        where(sc.eq.nan) sc=(/0.01,5.,34.3,0.,0./)
      case  (8)
        ! ssa, asym, hotpeak, hotwidth
        where(sc.eq.nan) sc=(/0.101,-0.263,0.589,0.046,0./) 
      case  (9)
        ! isotropic, volumetric, geometric, aspect1, aspect2
        where(sc.eq.nan) sc=(/1.,0.,0.,2.,1./)
      case default
        ! fraction of snow, ocean, sand, vegetation
        where(sc.eq.nan) sc=(/1.,0.,0.,0.,0./)
      end select

      
c assign user radiance angles
      
      if(radcalc) call vuangles(nphi,phi,nzen,uzen,nvzen,vzen,iout)

c check input parameters

      call chkin

c solar spectrum, atmosphric profile, satellite sensor
c filtering, surface albedo and aerosol parameters

      if(iday.ne.0) then
        call zensun(abs(iday),time,alat,alon,sza,saza,solfac)
      elseif(csza.ne.zip) then 
        sza=acos(csza)/dtor
      endif
      if(abs(sza-90).lt..01) sza=95.

      phi0=mod(saza-180.0_kr+360.0_kr, 360.0_kr)

      if(iday.lt.0) then 
        print '(a5,6a9)','day','time','lat','lon','sza','azm','solfac'
        print '(i5,6f9.3)',abs(iday),time,alat,alon,sza,saza,solfac
        if(radcalc) then
          print '(2a9)','phi','rel_az'
          do i=1,nphi 
            if(phi0.gt.180..and.phi(nphi)+phi0.gt.360) then
              print '(2f9.3)',phi(i),phi(i)+phi0-360.
            else
              print '(2f9.3)',phi(i),phi(i)+phi0
            endif
          enddo
        endif
        stop
      endif

      call setfilt(isat,wlinf,wlsup,wlinc,wl1,wl2,nwl)
      
      if (iout == 2) kdist = 0
      if(kdist.ge.0) then
        call atms(idatm,amix,nz,z,p,t,wh,wo)
        
        if(ngrid.ne.0.) call zgrid(nz,z,p,t,wh,wo,zgrid1,zgrid2,ngrid)

        if(zpres.ne.zip)  then
          call locate(z,nz,zpres,j)
          fj=(zpres-z(j))/(z(j+1)-z(j))
          pbar=p(j)*(p(j+1)/p(j))**fj
        endif

        call modatm(nz,sclh2o,uw,uo3,o3trp,ztrp,pbar,z,p,wh,wo)
        call modmix(xn2,xo2,xco2,xch4,xn2o,xco,xno2,xso2,xnh3,
     &       xno,xhno3,xo4)

        if(rhaer.lt.0.) rhaer=relhum(t(1),wh(1))
        
      elseif(kdist.eq.-1) then

        call gasinit(wl1,wl2,nwl,nz,z,p,t,rhaer,idb(3))

      endif
      
      if(maxval(idb).eq.0) call stdout0(iout,nwl,nz)


      temper(0)=t(nz)
      do j=1,nz
        temper(j)=t(nz+1-j)
      enddo
      if(btemp.lt.0.) btemp=temper(nz)
      if(ttemp.lt.0.) ttemp=temper(0)

      if(spowder) then          ! add a sub-surface bottom layer
        if(nz.ge.mxly) then
          print *,'Error --- nz < mxly is required with spowder option'
          stop
        endif
        z(2:nz+1)=z(1:nz)   ;  z(1)=-1.
        p(2:nz+1)=p(1:nz)   ;  p(1)=1.1*p(2)
        t(2:nz+1)=t(1:nz)   ;  t(1)=btemp 
        wh(2:nz+1)=wh(1:nz) ;  wh(1)=0.
        wo(2:nz+1)=wo(1:nz) ;  wo(1)=0.
        nz=nz+1
      endif

      nstrsv=nstr

      mcldz=max(numset(zero,tcloud,ncldz),
     &          numset(zero,lwp,ncldz))

      call zlayer(nz,z,mcldz,zcloud,lcld)

      if(kdist.ge.0) then 
        if(rhcld.ge.0) then
          if(krhclr.eq.1) then 
            call satcloud(nz,ncldz,lcld,t,rhcld,wh,idb(2))
          else
            call saturate(nz,ncldz,lcld,z,t,rhcld,wh,idb(2))
          endif
        endif
        if(ngrid.lt.0.or.idatm.lt.0) call prnatm
        call absint(uu,nz,z,p,t,wh,wo,idb(3))
      endif

c find nearest computational levels to given output altitudes
      
      if(minval(zout).lt.0.) then
        zout(1)=abs(zout(1))
        zout(2)=abs(zout(2))
        zflag=.true.
      endif
      call nearest(z,nz,zout(1),nbot)
      call nearest(z,nz,zout(2),ntop)

      nbot=nz-nbot+2
      ntop=nz-ntop+2
      if(ntop.eq.2) ntop=1      

      if(zflag) then
        print '(a,3i5,2f12.3)','nz,ntop,nbot,z(ntop),z(nbot): ',
     &       nz,ntop,nbot,z(nz-nbot+2),z(nz-ntop+2)
        stop
      endif

c set up for radiance calculation
  
      if(radcalc) then
        usrang=.true.
        numu=nzen
        do j=1,numu
          umu(j)=min(one,max(cos(uzen(numu+1-j)*dtor),-one))
          if(umu(j).eq.0.) then
            if(j.eq.numu) then 
              umu(j)=-.0001
            else
              umu(j)=.0001
            endif
          endif
        enddo
        if(nphi.ne.0) then
          azmavg=.false.
        else
          nphi=1
          phi(1)=0.
        endif
      endif

      if(idb(4).ne.0) call chkprn
c mfl
c     print *, 'drt: isalb = ', isalb
      call suralb(isalb,albcon,sc)

c load extinction, absorption, and asymmetry parameters for boundary
c layer aerosols, i.e., either rural, urban, oceanic or tropospheric aerosols
c and performs interpolation over relative humidity
      
      amu0=cos(sza*dtor)

c beginning of wavelength looping

      do                        ! wl_loop

        if(kdist.lt.0) then
          call readk(nz,wl1,wl2,wl,wvnmlo,wvnmhi,
     &         ib,nb,nk,etirr,ewcoef,gwk,dtauk,idb(7))
          if(nk.eq.0) exit      ! end of file, stop looping
        else
          call wllimits(nwl,wlinc,wl1,wl2,wl,wvnmhi,wvnmlo)
          if(nwl.eq.0) exit     ! wl beyond wl2, stop looping
          call gasset(kdist,wl,uu,amu0,nz,z,nk,gwk,dtauk,dtaugc,idb(7))
          ewcoef=1.
        endif

        dwl=10000./wvnmlo-10000./wvnmhi

        if (iout.eq.2) then
           if (first) print *, "nwl", nwl
           call taugas(wl,uu,amu0,nz,z,dtauc,dtauk(:,1),1)
           first=.false.
           cycle
        endif
           

        if(nf.ne.-2) etirr=solirr(wl,nf)*dwl

        flxin=etirr*solfac
        
        
c flxin set to et solar spectra (W/m2) or a uniform constant 
        
        if(nf.eq.0) flxin=dwl
        if(sza.ge.90.) then
          flxin=0.
          amu0=1.
        endif
        
        ff=filter(wl)*ewcoef
        
        if(nothrm.lt.0) then 
          plank=wl.gt.2.
        else
          plank=nothrm.eq.0
        endif
        
        if(any(isalb.eq.(/7,8,9/))) then
          lamber=.false.
          if(idb(8).ne.0) rsfc=dref(wvnmlo,wvnmhi,cos(dtor*sza))
          if(idb(8).eq.1) then
            print '(2f11.4)',wl,rsfc
          elseif(idb(8).eq.2) then
            call bdrefchk(wvnmlo,wvnmhi,sza,nphi,phi,nzen,uzen,rsfc)
          endif
        else
          lamber=.true.
          if(any(isalb.eq.(/-7,-8,-9/))) then
            rsfc=dref(wvnmlo,wvnmhi,cos(dtor*sza))
          else
            rsfc=salbedo(wl)
          endif
          rsfc=max(zero,(min(rsfc,one)))
          if(idb(8).ne.0) print '(2f11.4)',wl,rsfc
        endif
        
c calculate cloud optical depth and scattering parameters
        
        if(radcalc.and.corint) then  
          nmom=maxmom 
        else
          nmom=min(nstr+2,nstrms) ! add two to allow for nstr dithering
        endif
        
        dtaus=0.
        wreal=0.
        pmom=0.
        if(mcldz.gt.0) then
          call taucloud(nz,ncldz,wl,lcld,lwp,tcloud,nre,
     &         dtauc,wcld,imomc,nmom,pmom)
        elseif(nre(1).eq.zero) then 
          call usrcloud(nz,wl,p,dtauc,wcld,imomc,nmom,pmom,idb(5))
        endif

        if(idb(5).gt.0) call taucdb 
        
c aerosol optical depth
        
        call tauaero(wl,nz,z,nosct,dtaua,waer,nmom,pmom,idb(6))
        if(maxval(idb(1:8)).gt.0) cycle
        
c rayleigh scattering

        call rayleigh(wl,z,p,t,nz,dtaur)
        if(xrsc.ne.one) dtaur(1:nz)=xrsc*dtaur(1:nz)
        
c normalize phase function moments
        
        call normom(nz,dtauc,wcld,dtaua,waer,dtaur,nmom,pmom)
        
c       loop through k-distribution terms, calculate total optical
c       thickness (gases, rayleigh, aerosols and clouds)
        
        do kd=1,nk              ! kd_loop
          
          call depthscl(kdist, kd, nk, ib, nz, wl, dtaur, dtaua,
     &         waer, dtauc, wcld, spowder, gwk, dtauk, dtaugc, wt,
     &         dtaus, wreal, idb(9))
          if(idb(9).gt.0) cycle
          
          do j=0,2              ! retry loop in case radiance is negative
            nstr=nstrsv+j*(3*j-5)
            if(nstr.lt.4) cycle
            if(nstr.gt.nstrms.or.ff.eq.0.) exit
            
            call disort(nz,dtaus,wreal,corint,nmom,pmom,temper,
     &           wvnmlo,wvnmhi,usrtau,ntau,utau,nstr,usrang,
     &           numu,umu,nphi,phi,ibcnd,flxin,amu0,phi0,fisot,
     &           lamber,rsfc,btemp,ttemp,temis,plank,onlyfl,accur,
     &           prnt,header,mxly,maxulv,maxumu,maxphi,maxmom,
     &           rfldir,rfldn,flup,dfdt,uavg,uur,albmed,trnmed)
            
            if(nstr.gt.0) exit  ! success
          enddo
          
          if(nstr.lt.0.or.nstr.gt.nstrms) then
            write(*,*)'Error --- NSTR dithering procedure failed'
            stop
          endif
          nstr=nstrsv
          
          call stdout1(nz,z,ntop,nbot,iout,wl,dwl,wt,rfldir,rfldn,
     &         flup,ff,nphi,nzen,phi,uzen,uur,ib,nb,kd,nk)
          
        enddo                   ! kd_loop
      enddo                     ! wl_loop
      
      call stdout2(nz,iout,wl1,wl2,nphi,nzen,phi,uzen,z,p)
        
      contains 

c.......................................................................
      subroutine chkin          ! internal procedure

c PURPOSE     check input for invalid values. does not change anything.
c             if input is suspect, issue warnings and continue 
c             if input is bad print error messages and stop
c
c input:  a bunch of the namelist input parameters, by use association
c output: none
c
      integer :: kill
      character (len=132) :: line

                                ! warnings
      

      if(iaer.eq.0 .and. (vis.ne.zip.or.tbaer.ne.zip))
     &     call errmsg(16,'CHKIN--IAER=0, though VIS or TBAER set')

      if(corint .and. onlyfl)
     &     call errmsg(17,'CHKIN--CORINT=t, but flux output selected')

                                ! fatal errors
      kill=0
      
      if(idatm.lt.-6.or.idatm.gt.6) then
        call ck(kill,'idatm','[-6,6]')
        print *,'idatm=',idatm
      endif

      if(wlinf.lt.0.199) then
        call ck(kill,'wlinf','[0.2,-]')
        print *,'wlinf=',wlinf
      endif

      if(isat.le.-2) then
        if(wlsup.lt.0..or.wlsup.ge.wlinf) then
          call ck(kill,'wlsup','[0,wlinf]')
          print *,'wlsup=',wlsup
        endif
      else
        if(wlsup.lt.wlinf.or.wlsup.gt.100.) then
          call ck(kill,'wlsup','[wlinf,100]')
          print *,'wlsup=',wlsup
        endif
      endif

      if(isat.lt.-4.or.isat.gt.29) then
        call ck(kill,'isat','[-4,29]')
        print *,'isat=',isat
      endif

      if(solfac.lt.0.) then
        call ck(kill,'solfac','[0,inf]')
        print *,'solfac=',solfac
      endif

      if(minval(zcloud).lt.-100..or.maxval(zcloud).gt.100) then
        call ck(kill,'zcloud','[-100,100]')
        print *,'zcloud=',zcloud(1:ncldz)
      endif

      if(minval(abs(nre)).lt.2 .or. maxval(abs(nre)).gt.128) then
        if(nre(1).ne.0.) then
          call ck(kill,'nre','[2,128]')
          print *,'nre',nre(1:ncldz)
        endif
      endif

      if(any(tcloud.eq.0.and.zcloud.lt.0)) then
        print *,'CHKIN --- Error detected in input'
        print *,'TCLOUD(k)=0 when ZCLOUD(k)<0'
        print '(a,100f10.2)','zcloud',zcloud(1:ncldz)
        print '(a,100f10.2)','tcloud',tcloud(1:ncldz)
        kill=1
      endif

      if(minval(lwp).lt.0.) then
        call ck(kill,'lwp','[0,inf]')
        print *,'lwp=',lwp
      endif

      if(maxval(abs(zaer)).gt.100.) then
        call ck(kill,'zaer','[-100,100]')
        print *,'zaer=',zaer
      endif

      if(minval(taerst).lt.0.) then
        call ck(kill,'taerst','[0,inf]')
        print *,'taerst',taerst
      endif
        
      if(minval(jaer).lt.0.or.maxval(jaer).gt.4) then
        call ck(kill,'jaer','[0,4]')
        print *,'jaer',jaer
      endif
        
      if(nf.lt.-2.or.nf.gt.3) then
        call ck(kill,'nf','[-2,3]')
        print *,'nf',nf
      endif
        
      if(iaer.lt.-1.or.iaer.gt.5) then
        call ck(kill,'iaer','[-1,5]')
        print *,'iaer',iaer
      endif
        
      if(.not.any(isalb.eq.(/-7,-8,-9,-1,0,1,2,3,4,5,6,7,8,9,10/))) then
        call ck(kill,'isalb','[-7,-8,-9,-1,0,1,2,3,4,5,6,7,8,9,10]')
        print *,'isalb',isalb
      endif
        
      if(isalb.eq.0.and.albcon.lt.0.) then
        call ck(kill,'albcon','[0,inf]')
        print *,'albedo set by albcon',albcon
      endif

      if(minval(zout).lt.0..or.maxval(zout).gt.100)then
        call ck(kill,'zout','[0,100]')
        print *,'zout',zout
      endif
        
      if(.not. any(iout.eq.(/1,2,5,6,7,10,11,20,21,22,23/))) then
        call ck(kill,'iout','[1,2,5,6,7,10,11,20,21,22,23]')
        print *,'iout',iout
      endif
        
      if(nphi.lt.0.or.nphi.gt.nstrms) then
        call ck(kill,'nphi','[0,nstrms]')
        print *,'nphi',nphi
      endif
        

      if(any(iout.eq.(/20,21,22,23/)).and.nzen.eq.0) then
        write(*,'(1x,a,i2,a,a)') 'iout =',iout,
     &       ' implies radiance calculation, ',
     &       'but nzen=0 produces no radiance output'
        kill=1
      endif

      if(zpres.ne.zip.and.pbar.ne.zip) then
        write(*,'(1x,a)') 'set zpres or pbar but not both'
        write(*,*) 'zpres,pbar ',zpres,pbar
        kill=1
      endif

      if((numset(zero,tcloud,ncldz).ne.0).and.
     &   (numset(zero,lwp,ncldz).ne.0)) then
        write (*,*) 'set TCLOUD or LWP, but not both'
        kill=1
      endif

      if(kill.eq.1) then
        do 
          read(11,'(a)',end=100) line
          write(16,'(a)') trim(line)
        enddo
 100    continue
        stop
      endif
      return
      end subroutine chkin
c.......................................................................
      subroutine ck(kill,name,range) ! internal procedure (used in chkin)
      character*(*) name,range
      integer :: kill

      if(kill.eq.0) print '(a)','CHKIN --- Errors detected in INPUT'
      print '(/5x,4a)','Input parameter ',name,' not within ',range
      kill=1
      return
      end subroutine ck
c.......................................................................
      subroutine chkprn         ! internal procedure

c     purpose: print some info on variables defined within main program
c     does not change anything

      integer ::  i
      character*(30) :: form='("default value used for",a)'

      if(iday.ne.0) then
        print '(a,i3)',     'Day of year: ',abs(iday)
        print '(a,f8.4)',   '        GMT: ',time        
        print '(a,2f10.4)', '    lat,lon: ',alat,alon
        print *, ' '
      endif
      print '(a,2f10.4)', '  sza,solfac: ',sza,solfac
      print '(a,2f10.4)', '        zout: ',zout
      print '(a,3i10)',   'nbot,ntop,nz: ',nbot,ntop,nz
      print *, ' '
      if(pbar.lt.0)   print form,'PBAR'
      if(uw.lt.0)     print form,'UW'
      if(uo3.lt.0)    print form,'UO3'
      if(o3trp.lt.0)  print form,'O3TRP'
      if(sclh2o.lt.0) print form,'SCLH2O'
      print *, ' '
      if(tcloud(1).ne.0.) then
        call zlayer(nz,z,numset(zero,tcloud,ncldz),zcloud,lcld)
        print '(a,5f10.4)','    zcloud:' ,zcloud
        print '(a,5i10)',  '      lcld:' ,lcld
      endif
      
      print *,'"ATM'
      print *,nz
      print '(a/(10es11.3))', 'z',(z(i),i=1,nz)
      print '(a/(10es11.3))', 'p',(p(i),i=1,nz)
      print '(a/(10es11.3))', 't',(t(i),i=1,nz)
      print '(a/(10es11.3))', 'h2o',(wh(i),i=1,nz)
      print '(a/(10es11.3))', 'o3',(wo(i),i=1,nz)

      return

      end subroutine chkprn
c.......................................................................

      subroutine taucdb         ! internal procedure

! print cloud diagnostic

      integer :: ifirst=1

      if(ifirst.eq.1) print '(a,20i5)','lcld: ',lcld(1:mcldz)
      write(*,'(5x,4a11)',advance='no') 'wl','z','dtauc','ssa'
      do i=1,min(nmom,9)
        write(*,'("    pmom(",i1,")")',advance='no') i
      enddo
      print *
      do j=1,nz
        if(dtauc(j).eq.0) cycle
        print '(i5,13(f11.5))',j,wl,z(nz-j+1),dtauc(j),wcld(j),
     &       (pmom(i,j)/(dtauc(j)*wcld(j)),i=1,min(nmom,9))
      enddo

      end subroutine taucdb
c.......................................................................
      subroutine prnatm
      print *,nz
      do i=1,nz
        write(*,'(f11.3,4es11.3)') z(i),p(i),t(i),wh(i),wo(i)
      enddo
      stop
      end subroutine prnatm

      end program sbdart
c=======================================================================
      subroutine vuangles(nphi,phi,nzen,uzen,nvzen,vzen,iout) 
c
c     modifies values of nphi, phi, nzen, uzen
c
      use params, only: nstrms,  kr, zip

      real(kr) :: p1,p2,z1,z2,xxx,zipp
      real(kr) phi(nstrms),uzen(nstrms),vzen(nstrms)
      integer :: i,ii,nphi,nzen,nvzen,iout

      ! finds radiance viewing angles

      if(nphi.gt.0) then
        if(numset(zip,phi,nstrms).ne.2)
     &       write(*,'(a)') 'Error in MAIN -- ' //
     &       'must specify exactly 2 values of phi when nphi is set'
        if(nphi.gt.nstrms)
     &       write(*,'(a)') 'Error in Main -- ' //
     &       'specified nphi larger than nstrms'
        p1=min(phi(1),phi(2))
        p2=max(phi(1),phi(2))
        do i=1,nphi
          phi(i)=p1+(i-1)*(p2-p1)/float(nphi-1)
        enddo
      else
        nphi=numset(zip,phi,nstrms)
        if(nphi.eq.0) then      ! set default
          nphi=19
          p1=0.
          p2=180.
          do i=1,nphi
            phi(i)=p1+(i-1)*(p2-p1)/float(nphi-1)
          enddo
        endif
      endif
      
      zipp=90.                  ! convert to real(real_kind)
      nvzen=numset(zipp,vzen,nstrms)
      do i=1,nvzen
        uzen(i)=180.-vzen(i)
      enddo
      
      if(nzen.gt.0) then
        if(numset(zip,uzen,nstrms).ne.2)
     &       write(*,'(a)') 'Error in MAIN -- ' //
     &       'must specify exactly 2 values of uzen when nzen is set'
        if(nzen.gt.2*nstrms)
     &       write(*,'(a)') 'Error in Main -- ' //
     &       'specified nzen larger than nstrms'
        z1=min(uzen(1),uzen(2))
        z2=max(uzen(1),uzen(2))
        ii=0
        do i=1,nzen
          xxx=z1+(i-1)*(z2-z1)/float(nzen-1)
          if(abs(xxx-90.).gt..05) then
            ii=ii+1
            uzen(ii)=xxx
          endif
        enddo
        nzen=ii
      else
        nzen=numset(zip,uzen,nstrms)
        if(nzen.eq.0) then      ! set default radiance grid
          select case (iout)
          case (5,20)  ; nzen=18 ; z1=0.  ; z2=85. ! toa looking down
          case (6,21)  ; nzen=18 ; z1=95. ; z2=180. ! surface looking up
          case default ; nzen=36 ; z1=0.  ; z2=180.
          end select
          do i=1,nzen
            uzen(i)=z1+(z2-z1)*(i-1)/float(nzen-1)
          enddo
          ! assume a reasonable number of streams if nstr unset by user
          if(nstr.eq.4) nstr=min( 2*(max(nphi,nzen)/2), nstrms)
        endif
      endif
      return
      end subroutine vuangles
        
c=======================================================================
      subroutine stdout0(iout,nwl,nz)

      implicit none
      integer :: iout, nwl, nz

      select case (iout)
      case (1,5,6) 
        write(*,'(/,a)') '"tbf'
        write(*,'(i15)') nwl
      case (7) 
        write(*,'(/,a)') '"fzw' 
        write(*,'(i15)') nz
      end select

      return
      end
c=======================================================================
      subroutine stdout1(nz,z,ntop,nbot,iout,wl,dwl,wt,
     &           rfldir,rfldn,flup,ff,nphi,nzen,phi,uzen,uur,
     &           ib,nb,kd,nk)
c
c performs these functions (depending on value of iout):
c  input:
c    nz       number of levels
c    z        altitude array [z(1) at bottom]
c    ntop     output level top
c    nbot     output level bottom
c    iout
c          1   print output at each wavelength (fluxes in w/m2/micron)
c          5   print radiance output at top at each wavelength
c          6   print radiance output at bot at each wavelength
c         10   total inband flux (w/m2) summed over wavelength 
c         11   flux at each level integrated over wavelength
c         20   total inband radiance (w/m2/str) at ZOUT(2)
c         21   total inband radiance (w/m2/str) at ZOUT(1)
c         22   inband radiance (w/m2/str) at each level
c         23   inband radiance (w/m2/str) at zout(1) and zout(2)
c                upper hemisphere radiance at zout(2)
c                lower hemisphere radiance at zout(1)
c    wl      current wavelength (um)
c    dwl     wavelength increment (um)
c    wt      weighting factor for this iteration
c    rfldir  direct flux
c    rfldn   downward flux
c    flup    upward flux
c    ff      filter function
c    nphi    number of azimuth angles in radiance viewing array
c    nzen    number of zenith angles in radiance viewing array
c    phi     azimuth viewing angle array
c    uzen    zenith viewing angle array
c    uur     radiance array
c    ib      band index.
c    nb      total number of sub-bands
c    kd      k-distribution index
c    nk      number k-distribution intervals
c  
      use params, only: mxly, kr, nstrms
      use outblk, only: topdn,topup,topdir,botdn,botup,botdir,
     &     maxulv, uurs, phidw, fxdn, fxup, fxdir, uurl
      implicit none

      integer, parameter :: nr=mxly*3, nta=9

      integer :: iout, kd, ntop, nbot, nk, i, j, k, nphi, nzen,
     &     nz, im, ib, nb

      real(kr) :: rfldir(*),rfldn(*),flup(*),phi(*),uzen(*),z(*),
     &     uur(nstrms,maxulv,nstrms), dwt, wt, dwl, ff, wl, 
     &     dwx, weq, wfull

 1000 format(10es12.4)

      dwt=wt*ff       ! flxin in w/m2

      if(iout.eq.1 .or. iout.eq.5 .or. iout.eq.6) then
        if(kd.eq.1.and.ib.eq.nb) then
          topdn=0.
          topup=0.
          topdir=0.
          botdn=0.
          botup=0.
          botdir=0.
          weq=0. 
          wfull=0.
        endif
        topdn  = topdn  + (rfldn(ntop)+rfldir(ntop))*dwt
        topup  = topup  + flup(ntop)*dwt
        topdir = topdir + rfldir(ntop)*dwt
        botdn  = botdn  + (rfldn(nbot)+rfldir(nbot))*dwt
        botup  = botup  + flup(nbot)*dwt
        botdir = botdir + rfldir(nbot)*dwt

        if(kd.eq.nk) then
          weq=weq+dwl*ff        ! finished sub-band, increment equiv width
          wfull=wfull+dwl       ! finished sub-band, increment full width
        endif

        if(kd.eq.nk.and.ib.eq.1) then
          if(weq.eq.0.) weq=1.e-30
          write(*,'(f12.8,f9.5,6es12.4)') wl,weq/wfull,
     &       real(topdn/weq),real(topup/weq),real(topdir/weq),
     &       real(botdn/weq),real(botup/weq),real(botdir/weq)
        endif

        if(iout.eq.5 .or. iout.eq.6 ) then
          j=ntop
          if(iout.eq.6) j=nbot

          if(kd.eq.1.and.ib.eq.nb) uurs(1:nzen,1:nphi)=0.
          
          do  k=1,nphi
            do  i=1,nzen
              uurs(i,k)=uurs(i,k)+uur(i,j,k)*dwt
            enddo
          enddo

          if(kd.eq.nk.and.ib.eq.1) then
            if(weq.eq.0.) weq=1.e-30
            write(*,'(3i4)')       nphi,nzen
            write(*,1000) (real(phi(j)),j=1,nphi)
            write(*,1000) (real(uzen(j)),j=1,nzen)
            do i=nzen,1,-1
              write(*,1000) (real(uurs(i,k)/weq),k=1,nphi)
            enddo
          endif
        endif
      endif

      if(any(iout.eq.(/10,11,20,21,22,23/)) .and. kd.eq.nk)
     &     phidw=phidw+dwl*ff

      if(iout.eq.7 .or. iout.eq.11 .or. iout.eq.22) then
        if(iout.eq.7) then
          if(kd.eq.1.and.ib.eq.nb) then
            fxdn(1:nz)=0.
            fxup(1:nz)=0.
            fxdir(1:nz)=0.
          endif
        endif
        do 10 i=1,nz
          fxdn(i)=fxdn(i)+(rfldn(i+1)+rfldir(i+1))*dwt
          fxup(i)=fxup(i)+flup(i+1)*dwt
          fxdir(i)=fxdir(i)+rfldir(i+1)*dwt
 10     continue
      endif

      if(iout.eq.7.and.kd.eq.nk.and.ib.eq.1) then
        write(*,'(//,f12.8)') wl
        write(*,'(/(10es11.3))') (z(im),im=nz,1,-1)
        write(*,'(/(10es11.3))') (real(fxdir(i)),i=1,nz)
        write(*,'(/(10es11.3))') (real(fxdn(i)-fxdir(i)),i=1,nz)
        write(*,'(/(10es11.3))') (real(fxdn(i)),i=1,nz)
        write(*,'(/(10es11.3))') (real(fxup(i)),i=1,nz)
      endif

      if(any(iout.eq.(/10,20,21,23/))) then
        topdn  = topdn  + (rfldn(ntop)+rfldir(ntop))*dwt
        topup  = topup  + flup(ntop)*dwt
        topdir = topdir + rfldir(ntop)*dwt
        botdn  = botdn  + (rfldn(nbot)+rfldir(nbot))*dwt
        botup  = botup  + flup(nbot)*dwt
        botdir = botdir + rfldir(nbot)*dwt
      endif

      if(iout.eq.20 .or. iout.eq.21) then
        j=ntop
        if(iout.eq.21) j=nbot
        do 30 k=1,nphi
          do 20 i=1,nzen
            uurs(i,k)=uurs(i,k)+uur(i,j,k)*dwt
 20       continue
 30     continue
      endif

      if(iout.eq.22) then
        do 60 j=1,nz
          do 50 k=1,nphi
            do 40 i=1,nzen
              uurl(k,i,j)=uurl(k,i,j)+uur(i,j+1,k)*dwt
 40         continue
 50       continue
 60     continue
      endif
      
      if(iout.eq.23) then
        do 80 k=1,nphi
          do 70 i=1,nzen
            if(uzen(nzen-i+1).lt.90.) then           
              j=ntop
            else
              j=nbot
            endif
            uurs(i,k)=uurs(i,k)+uur(i,j,k)*dwt
 70       continue
 80     continue
      endif

      return

      end

c======================================================================
      subroutine stdout2(nz,iout,wlinf,wlsup,nphi,nzen,phi,uzen,
     &     z,p)

      use params, only: grav, kr
      use outblk, only: phidw, fxdn, fxup, fxdir, topdn, topup, topdir,
     &     botdn, botup, botdir, uurs, uurl
      implicit none

      integer, parameter :: nta=9

      real(kr), parameter :: cp=1004. ! heat cap of air at constant p (J/K/kg)
      
      integer :: iout, j, nz, i, nphi, nzen, k

      real(kr) :: z(*), p(*), phi(*),
     &     uzen(*), fntm, zz, pp, fnt, dfdz, heat, wlinf, wlsup,
     &     zm, pm, eps=1.e-30

      character*(72) :: outstr
c
 1000 format(10es12.4)

      if(iout.eq.11) then
        write(*,'(i4,es15.7)') nz,phidw
c       fxdn(nz) corresponds to total downward flux at surface
        fntm=0.
        do 10 i=1,nz
          zz=z(nz-i+1)
          pp=p(nz-i+1)
          fnt=fxdn(i)-fxup(i)
          if(i.gt.1) then
            dfdz=(fntm-fnt)/(zm-zz)
            heat=.01*3600*24*grav*(fntm-fnt)/(cp*(pp-pm)) ! kelvin/day
          else
            dfdz=0.
            heat=0.
          endif
          fntm=fnt
          zm=zz
          pm=pp
          write(*,1000) zz,pp,real(fxdn(i)),real(fxup(i)),
     &         real(fxdir(i)),real(dfdz),real(heat)
 10     continue
      endif

      if(iout.eq.10.or.iout.eq.20.or.iout.eq.21.or.iout.eq.23)
     &   write(*,'(3f11.4,6es12.4)') wlinf,wlsup,phidw,
     &     real(topdn),real(topup),real(topdir),
     &     real(botdn),real(botup),real(botdir)

      if(iout.eq.20.or.iout.eq.21.or.iout.eq.23) then
        write(*,'(3i4)')       nphi,nzen
        write(*,1000) (phi(j),j=1,nphi)
        write(*,1000) (uzen(j),j=1,nzen)
        do 20 i=nzen,1,-1
          write(*,'(20es12.4)' ) (real(uurs(i,k)),k=1,nphi)
 20     continue
      endif

      if(iout.eq.22) then
        write(*,'(3i4,es12.4)') nphi,nzen,nz,phidw
        write(*,1000) (phi(i),i=1,nphi)
        write(*,1000) (uzen(j),j=1,nzen)
        write(*,1000) (z(k),k=nz,1,-1)
        write(*,1000) (real(fxdn(k)),k=1,nz)
        write(*,1000) (real(fxup(k)),k=1,nz)
        write(*,1000) (real(fxdir(k)),k=1,nz)
        write(*,1000)
     &       (((real(uurl(i,j,k)),i=1,nphi),j=nzen,1,-1),k=1,nz)
      endif
      return
      end
c=======================================================================
c purpose: returns the largest index i, such zz(i).ne.zip

      integer function numset(zip,zz,nz)
      use params, only: kr
      implicit none
      integer :: nz,i
      real(kr) :: zz(nz),zip
      numset=0
      do 10 i=1,nz
        if(zz(i).ne.zip) numset=i
 10   continue
      return
      end      
c=======================================================================
      subroutine locate(xx,n,x,j)
c
c purpose:  given an array xx of length n, and given a value X, returns
c           a value J such that X is between xx(j) and xx(j+1). xx must
c           be monotonic, either increasing of decreasing. this function
c           returns j=1 or j=n-1 if x is out of range. 
c
c input:
c   xx      monitonic table
c   n       size of xx
c   x       single floating point value perhaps within the range of xx
c
c output:
c           function returns index value j, such that 
c          
c            for an increasing table
c
c                xx(j) .lt. x .le. xx(j+1),  
c                j=1    for x .lt. xx(1)
c                j=n-1  for x .gt. xx(n)
c
c            for a decreasing table
c                xx(j) .le. x .lt. xx(j+1)
c                j=n-1  for x .lt. xx(n)
c                j=1    for x .gt. xx(1)
c
      use params, only: kr
      implicit none
      integer :: j, n, jl, jm, ju
      real(kr) :: x,xx(n)

c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      if(x.eq.xx(1)) then 
        j=1
        return
      endif
      if(x.eq.xx(n)) then
        j=n-1
        return
      endif
      jl=1
      ju=n
10    if(ju-jl.gt.1) then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      end
c=======================================================================
      subroutine nearest(xx,n,x,j)
c     purpose: for array xx of size n, find j such that xx(j) is nearest x

      use params, only: kr
      implicit none
      integer :: j, k, n
      real(kr) :: xx(*), x
      call locate(xx,n,x,k)
      j=k
      if(abs(x-xx(k+1)).lt.abs(x-xx(k))) j=k+1
      return
      end
c=======================================================================
c purpose: for array xx of size nz, return xx in reverse order

      subroutine reverse(nz,xx)

      use params, only: kr
      implicit none
      integer :: i, nz, ir
      real(kr) :: xx(*), x
      do 10 i=1,nz/2
        ir=nz-i+1
        x=xx(i)
        xx(i)=xx(ir)
        xx(ir)=x
 10   continue
      return
      end
c=======================================================================

      subroutine zlayer(nz,z,nzz,zz,lz)
c   
c purpose: find the index array lz(k) such that z(nz-lz(k)+1) is the
c          largest element of z less than or equal to zz(i)+.001
c          (note reversal of array order).  if zz(k) < 0 for k>1 then
c          lz(k) will also be tagged negative. Other routines use this info
c          to indicate altitude ranges.  the smallest value of z
c          is normally zero, but the first element of z can be negative
c          to indicate one subsurface layer.
c
c          REMEMBER: z, t, p, wh, wo all use bottom up numbering while
c                    everthing else uses top-down numbering 
c

c input:
c   nz     number of atmospheric layers
c
c   z      array of layer altitudes (monotonically increasing)
c
c   nzz    number of flagged altitude layers
c
c   zz     array of flagged altitudes
c
c
c output:
c lz       index array of flagged altitudes 
c
c example: nz = 9 
c                   9    8   7   6    5   4   3   2   1
c          z  = [   0,   1,  3,  5,   7,  9, 10, 12, 15 ]
c          zz = [   0,   2, -3,  3,  -8]
c          lz = [   9,   8,   8,   7, 5]
c
      use params, only: kr

      implicit none
      integer :: lz(*), kk, isgn, j, jj, nz, nzz
      real(kr) :: z(*),zz(*),zcmpr
c

      do kk=1,nzz

        isgn=1
        if(kk.eq.1) then
          zcmpr=zz(kk)+.001
        else
          if(zz(kk).lt.0.) isgn=-1
          zcmpr=abs(zz(kk)+.001)
        endif
        
        do j=nz,1,-1
          if(z(j).le.zcmpr) goto 10
        enddo
        isgn=0
 10     continue
        lz(kk)=isgn*(nz-j+1)
      enddo

      return
      end
c=======================================================================
      subroutine levrng(nzz,lz,n,lbot,ltop)
c
c     purpose: used to specify a computational level range on the 
c              the flag array lz
c     
      implicit none
      integer :: lz(*), lbot, ltop, nzz, n

      lbot=0
      ltop=0
      if(n.gt.nzz) return
      if(lz(n).le.0) return

      lbot=lz(n)
      ltop=lbot
      if(n.ne.nzz.and.lz(n+1).lt.0) ltop=-lz(n+1)

      return
      end
c===========================================================================
      function planck(wl,t)

c     compute planck function at wavelength wl and temperature t.  
c     units are w/m2/sr/micron (note that w/m2/sr/micron = w/cm2/sr/cm)

      use params, only: kr
      implicit none

      real(kr), parameter :: h=6.6262e-27, bk=1.3807e-16, c=2.9979e10,
     &                   c1=1.e13*2*h*c**2, c2=1.e4*h*c/bk
      real(kr) :: wl, t, planck
c 
      planck=0.
      if(wl.gt.1.) planck=c1/(wl**5*(exp(c2/(wl*t))-1.))
      return
      end
c=======================================================================
      subroutine normom(nz,dtauc,wcld,dtaua,waer,dtaur,nmom,pmom)
c
c input:
c   nz      number of computational levels
c   dtauc   optical depth increment due to clouds
c   wcld    single scattering albedo of cloud
c   dtaua   optical depth increment due to aerosols
c   waer    single scattering albedo of cloud
c   dtaur   optical depth increment due to rayleigh scattering
c   nmom    number of legendre terms in pmom
c   pmom    phase function moments weighted by scattering
c           optical depth (clouds and aerosol only)
c
c output:
c   pmom    normalized phase function moments due to clouds, 
c           aerosol and rayleigh
c
      use params, only: maxmom,kr
      implicit none

      integer :: i,nz,nmom
      real(kr) :: dtauc(*),wcld(*),dtaua(*),waer(*),
     &        dtaur(*),pmom(0:maxmom,*), dtsct

      do i=1,nz
        pmom(2,i)=pmom(2,i)+.1*dtaur(i)
        dtsct=dtauc(i)*wcld(i)+dtaua(i)*waer(i)+dtaur(i)
        if(dtsct.ne.0.) pmom(0:nmom,i)=pmom(0:nmom,i)/dtsct
      enddo
      pmom(0,1:nz)=1.
      return
      end
c=======================================================================
      subroutine bdrefchk(wv1,wv2,sza,nphi,phi,nzen,uzen,albedo)

c     print the bdrf function in the iout=20 format.  This diagnostic 
c     output may be compared to the radiance profile produced with inputs
c     (iout=20, pbar=0, isalb=7,8,9),  which will yield the 
c     radiance from the same bdrf function fit by legendre polynomials.  
c     The bdrf comparison should make use of the defining equation of bdrf:

c        bdrf=pi*R/botdn 

c     another diagnostic output produced by this routine is the 
c     albedo produces by the test bdrf.  This quatity is listed in the two 
c     iout=20 output positions usually used for botup and topup.  The 
c     slots filled by topdn, topdir, botdn, and botdir are set to one.

      use params, only: pi,kr,one
      implicit none
      integer :: i,j,k,nphi,nzen,idb
      real(kr) :: dtor, amu0, bdref, zenr, sza, phi(nphi),uzen(nzen),
     &     amu,dmu,albedo,bdrf(nphi,nzen),work(nzen),wl,wv1,wv2
      character (len=8) :: fmt='(..f...)'

      wl=20000./(wv1+wv2)
      dtor=pi/180.
      amu0=cos(dtor*(sza))
      do j=1,nzen
        do i=1,nphi
          bdrf(i,j)=bdref(wv1,wv2,amu0,cos(dtor*uzen(j)),dtor*phi(i))
        enddo
      enddo

      print '(9f8.3)',wl,0.,0.,0.,0.,0.,1.,albedo,0.
      print '(2i5)',nphi,nzen
      
      fmt(5:7)='8.2'
      write(fmt(2:3),'(i2)') nphi ; print fmt,phi(1:nphi)
      write(fmt(2:3),'(i2)') nzen ; print fmt,uzen(1:nzen)
      fmt(5:7)='8.4'
      write(fmt(2:3),'(i2)') nphi
      print fmt,bdrf(1:nphi,1:nzen)

      return
      end
c=======================================================================
      subroutine helper(iout)
      implicit none
      integer :: iout

      select case (iout)
      case (1) 
        print '(a)',
     &   'iout=1,  one output record for each wavelength, ',
     &   'output quantities are,',
     &   '',
     &   'WL,FFV,TOPDN,TOPUP,TOPDIR,BOTDN,BOTUP,BOTDIR',
     &   '',
     &   'WL    = wavelength                         (microns)',
     &   'FFV   = filter function value',
     &   'TOPDN = total downward flux at ZOUT(2) km  (w/m2/micron)',
     &   'TOPUP = total upward flux at ZOUT(2) km    (w/m2/micron)',
     &   'TOPDIR= direct downward flux at ZOUT(2) km (w/m2/micron)',
     &   'BOTDN = total downward flux at ZOUT(1) km  (w/m2/micron)',
     &   'BOTUP = total upward flux at  ZOUT(1) km   (w/m2/micron)',
     &   'BOTDIR= direct downward flux at ZOUT(1) km (w/m2/micron)',
     &   '',
     &   'NOTE: When ISAT ne 1 these radiometric quantities',
     &   '      are each multiplied by the filter function',
     &   '      To get the actual specific irradiance divide',
     &   '      by FFV(WL).',
     &   ''
        
      case (5,6)
        print '(a)',
     &   'iout=5 (output for zout(2) or iout=6, output for zout(1)',
     &   '(nzen+3) records for each wavelength.  Output format:',
     &   ' write(*,*) ''"tbf''  ; Block id (used in postprocessors)',
     &   '',
     &   '   do m=1,nw',
     &   '     write(*,*)',
     &   '  &    wl,ffv,topdn,topup,topdir,botdn,botup,botdir',
     &   '     write(*,*)  nphi,nzen',
     &   '     write(*,*) (phi(j),j=1,nphi)',
     &   '     write(*,*) (uzen(j),j=1,nzen)',
     &   '     do i=nzen,1,-1',
     &   '       write(*,*) (uurs(i,k),k=1,nphi)',
     &   '     enddo',
     &   '   enddo',
     &   '',
     &   ' where,',
     &   '',
     &   ' WL    = wavelength                         (microns)',
     &   ' FFV   = filter function value',
     &   ' TOPDN = total downward flux at ZOUT(2) km  (w/m2/micron)',
     &   ' TOPUP = total upward flux at ZOUT(2) km    (w/m2/micron)',
     &   ' TOPDIR= direct downward flux at ZOUT(2) km (w/m2/micron)',
     &   ' BOTDN = total downward flux at ZOUT(1) km  (w/m2/micron)',
     &   ' BOTUP = total upward flux at  ZOUT(1) km   (w/m2/micron)',
     &   ' BOTDIR= direct downward flux at ZOUT(1) km (w/m2/micron)',
     &   ' NPHI  = number of user azimuth angles',
     &   ' NZEN  = number of user zenith angles',
     &   ' PHI   = user specified azimuth angles      (degrees)',
     &   ' UZEN  = user specified zenith angles       (degrees)',
     &   ' UURS  = radiance at user angles at         (w/m2/um/str) ',
     &   '         altitude ZOUT(2) (top)',
     &   '',
     &   '',
     &   ' NOTE: The radiance output from SBDART represents',
     &   '       scattered radiation. It does not include the',
     &   '       solar direct beam.  Also, keep in mind that UURS',
     &   '       represents the radiance at the user specified',
     &   '       sample directions.  Hence, computing the',
     &   '       irradiance by an angular integration of UURS',
     &   '       will not yield BOTDN because of the neglect of',
     &   '       the direct beam, and it will probably not yield',
     &   '       (BOTDN-BOTDIR) because of under-sampling.',
     &   ''
        
      case (7)
        print '(a)','iout=7',
     &   'radiative flux at each layer for each wavelength.',
     &   '',
     &   'write(*,*) ''"fzw''     ; block id (used in postprocessors)',
     &   'write(*,*) nz           ; number of z levels',
     &   'write(*,*) nw           ; number of wavelengths',
     &   '',
     &   ' do j=1,nw',
     &   '   write(*,*) wl',
     &   '   write(*,*) ',
     &   '&    (Z(i),i=nz,1,-1),   ; altitude              (km)',
     &   '&    (fdird(i),i=1,nz),  ; downward direct flux  (w/m2/um)',
     &   '&    (fdifd(i),i=1,nz),  ; downward diffuse flux (w/m2/um)',
     &   '&    (flxdn(i),i=1,nz),  ; total downward flux   (w/m2/um)',
     &   '&    (flxup(i),i=1,nz)   ; total upward flux     (w/m2/um)',
     &   ' enddo',
     &   ''
        
      case (10)
        print '(a)','iout=10',
     &   'one output record per run, integrated over wavelength.',
     &   'output quantities are, (integrations by trapezoid rule)',
     &   '',
     &   '   WLINF,WLSUP,FFEW,TOPDN,TOPUP,TOPDIR,BOTDN,BOTUP,BOTDIR',
     &   '',
     &   '   WLINF = lower wavelength limit             (microns)',
     &   '   WLSUP = upper wavelength limit             (microns)',
     &   '   FFEW  = filter function equivalent width   (microns)',
     &   '   TOPDN = total downward flux at ZOUT(2) km  (w/m2)',
     &   '   TOPUP = total upward flux at ZOUT(2) km    (w/m2)',
     &   '   TOPDIR= direct downward flux at ZOUT(2) km (w/m2)',
     &   '   BOTDN = total downward flux at ZOUT(1) km  (w/m2)',
     &   '   BOTUP = total upward flux at  ZOUT(1) km   (w/m2)',
     &   '   BOTDIR= direct downward flux at ZOUT(1) km (w/m2)',
     &   ''
        
      case (11)
        print '(a)','iout=11',
     &   'radiant fluxes at each atmospheric layer integrated ',
     &   'over wavelength. Output format:',
     &   '',
     &   '   write(*,*) nz,phidw',
     &   '   do i=1,nz',
     &   '     write(*,*) zz,pp,fxdn(i),fxup(i),fxdir(i),dfdz,heat',
     &   '   enddo',
     &   '',
     &   ' where,  nz    = number of atmospheric layers',
     &   '         ffew  = filter function equivalent width(um)',
     &   '         zz    = level altitudes                 (km)',
     &   '         pp    = level pressure                  (mb)',
     &   '         fxdn  = downward flux (direct+diffuse)  (W/m2)',
     &   '         fxup  = upward flux                     (W/m2)',
     &   '         fxdir = downward flux, direct beam only (W/m2)',
     &   '         dfdz  = radiant energy flux divergence  (mW/m3)',
     &   '         heat  = heating rate                    (K/day)',
     &   '',
     &   'NOTE:  dfdz(i) and heat(i) are defined at the layer ',
     &   'centers, i.e., halfway between level i-1 and level i.',
     &   ''
        
      case (20)
        print '(a)',
     &   'iout=20,  radiance output at ZOUT(2) km. ',
     &   '',
     &   'Output format:',
     &   '',
     &   '   write(*,*)  wlinf,wlsup,ffew,topdn,topup,topdir,',
     &   '  &        botdn,botup,botdir',
     &   '   write(*,*)  nphi,nzen',
     &   '   write(*,*)  (phi(i),i=1,nphi)',
     &   '   write(*,*)  (uzen(j),j=1,nzen)',
     &   '   write(*,*)  ((r(i,j),i=1,nphi),j=1,nzen)',
     &   '',
     &   'The first record of output is the same as format IOUT=10',
     &   '(WLINF,WLSUP,FFEW,TOPDN,TOPUP,TOPDIR,BOTDN,BOTUP,BOTDIR)',
     &   'addition records contain:',
     &   '',
     &   '   NPHI  = number of user azimuth angles',
     &   '   NZEN  = number of user zenith angles',
     &   '   PHI   = user azimuth angles (nphi values) ',
     &   '   UZEN  = user zenith angles (nzen values)',
     &   '   R     = radiance array (nphi,nzen) (W/m2/sr)',
     &   ''
        
      case (21)
        print '(a)',
     &   'iout=21,  radiance output at ZOUT(1) km.',
     &   '',
     &   'Output format:',
     &   '',
     &   '   write(*,*)  wlinf,wlsup,ffew,topdn,topup,topdir,',
     &   '  &        botdn,botup,botdir',
     &   '   write(*,*)  nphi,nzen',
     &   '   write(*,*)  (phi(i),i=1,nphi)',
     &   '   write(*,*)  (uzen(j),j=1,nzen)',
     &   '   write(*,*)  ((r(i,j),i=1,nphi),j=1,nzen)',
     &   '',
     &   'The first record of output is the same as format IOUT=10',
     &   '(WLINF,WLSUP,FFEW,TOPDN,TOPUP,TOPDIR,BOTDN,BOTUP,BOTDIR)',
     &   'addition records contain:',
     &   '',
     &   '   NPHI  = number of user azimuth angles',
     &   '   NZEN  = number of user zenith angles',
     &   '   PHI   = user azimuth angles (nphi values) ',
     &   '   UZEN  = user zenith angles (nzen values)',
     &   '   R     = radiance array (nphi,nzen) (W/m2/sr)',
     &   ''
        
      case (22)
        print '(a)',
     &   'iout=22,  radiance and flux at each atmospheric layer ',
     &   'integrated over wavelength. ',
     &   '',
     &   'Output format:',
     &   '',
     &   '   write(*,*) nphi,nzen,nz,ffew',
     &   '   write(*,*) (phi(i),i=1,nphi)',
     &   '   write(*,*) (uzen(j),j=1,nzen)',
     &   '   write(*,*) (z(k),k=nz,1,-1)',
     &   '   write(*,*) (fxdn(k),k=1,nz)',
     &   '   write(*,*) (fxup(k),k=1,nz)',
     &   '   write(*,*) (fxdir(k),k=1,nz)',
     &   '   write(*,*) (((uurl(i,j,k),i=1,nphi),j=1,nzen),k=1,nz)',
     &   '',
     &   'where,  nphi  = number of user specified azimuth angles',
     &   '        nzen  = number of user specified zenith angles',
     &   '        nz    = number of atmospheric levels',
     &   '        ffew  = filter function equivalent width (um)',
     &   '        phi   = user specified anizmuth angles   (degrees)',
     &   '        uzen  = user specified zenith angles     (degrees)',
     &   '        z     = altitudes of atmospheric layers  (km)',
     &   '        fxdn  = downward flux (direct+diffuse)   (W/m2)',
     &   '        fxup  = upward flux                      (W/m2)',
     &   '        fxdir = downward flux, direct beam only  (W/m2)',
     &   '        uurl  = radiance at each layer           (W/m2/str)',
     &   ''
      end select
      return
      end
c=======================================================================
      subroutine wllimits(nwl,wlinc,wl1,wl2,wl,wvnmhi,wvnmlo)

c     purpose: compute wavelength for a given iteration index
c     input:
c       nwl    number of wavelength steps
c       wlinc  wavelength increment (used when kdist ge 0)
c
c     output:
c       nwl    set to zero when il exceeds nwl
c       wl     wavelength (um)
c       wvnmhi wavenumber upper limit used in disort (cm-1)
c       wnvmlo wavenumber lower limit used in disort (cm-1)
c
c     saved:
c       il     wavelength/wavenumber iteration index (each call to
c              wllimits causes il to be incremented by one)

      use params, only: kr

      implicit none
      integer :: il=0, nwl, kdist
      real(kr) :: wlinc,wl,wvnmhi,wvnmlo,dw,dwl,wi,wl1,wl2
      real(kr) :: ww1,ww2,wr,wvnm,dwvn

      if(il.ge.nwl) then
        nwl=0
        return
      endif

      wi=float(il)
      if(wlinc.gt.1) then       ! equal increments of wavenumber
        wl=wlwn(wi)
        ww1=wlwn(wi-.5)
        ww2=wlwn(wi+.5)
      elseif(wlinc.lt.0.) then  ! equal increments of log wavelength
        wr=wl2/wl1
        wl=wlln(wi)
        ww1=wlln(wi-.5)
        ww2=wlln(wi+.5)
      else                      ! equal increments of wavelength
        wl=wl1+wi*wlinc
        ww1=wl-.5*wlinc
        ww2=wl+.5*wlinc
      endif
      
      if(il.eq.0.and.il.ne.nwl-1)     ww1=wl
      if(il.eq.nwl-1.and.il.ne.0)     ww2=wl

      if(ww1.eq.wl .and. ww2.eq.wl) then
        ww1=wl-.0005
        ww2=wl+.0005
      endif

      wvnmlo=10000./ww2
      wvnmhi=10000./ww1

c      if(il.eq.0) print '(a5,7a11)',
c     &     'il','wlinc','ww1','ww2','wl','wvlo','wvhi','dw'
c      print '(i5,7es11.3)', il,wlinc,ww1,ww2,wl,wvnmlo,wvnmhi,
c     &     10000./wvnmlo-10000./wvnmhi

      il=il+1

      return      
      contains
c.......................................................................
      real(kr) function wlwn(x)          ! internal procedure
                                ! equal increments of wavenumber
      real(kr) :: x, xx

      xx=x/(nwl-1)
      wlwn=wl1*wl2/((1.-xx)*wl2+xx*wl1)
      return
      end function wlwn
c.......................................................................
      real(kr) function wlln(x)          ! internal procedure
                                ! equal increments of log wavelength
      real(kr) :: x, xx
      xx=x/(nwl-1)
      wlln=wl1*wr**xx
      return
      end function wlln

      end subroutine wllimits

