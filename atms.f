c file:                  atms.f
c
c external routines:     atms,modatm,saturate,satcloud
c
c internal routines      midsum,midwin,subsum,subwin,tropic,us62,useratm
c                        zgrid
c
c logical units:         unit 13, used to read user atmosphere from atms.dat

c=======================================================================
      subroutine satcloud(nz,ncldz,lcld,t,rhcld,wh,idb)
      use params, only: kr, tzero
      implicit none
      integer :: j, jj, lbot, ltop, idb,i,nz,ncldz
      real(kr) ::  t(*), wh(*), lcld(*), a, rhcld
c
c purpose: modify the watervapor density inside clouds to have a
c          relative humidity of RHCLD.
c
c INPUT:
c   nz        number of atmospheric levels
c   ncldz     size cloud layer arrays
c   lcld      cloud layer array
c   t         temperature at altitude z
c   rhcld     relative humidity of cloud layer, used to specify
c             water vapor density.
c
c
c h2osat is the mass density (g/m3) of water vapor at 100% saturation
c (source: handbook of chemistry and physics, h2o vapor pressure
c table d-112, assuming density related to pressure by ideal gas law)
c

      if(idb.gt.0) print '(a/(10es11.3))', 'wh_in  ',
     &           (wh(i)/satden(tzero/t(i)),i=1,nz)
c
c  adjust cloud layer, don't adjust clear layers
c
      do i=1,ncldz
        call levrng(ncldz,lcld,i,lbot,ltop)
        if(lbot.ne.0) then
          do j=max(ltop-1,1),lbot ! saturate to top of cloud
            jj=nz-j+1
            wh(jj)=rhcld*satden(tzero/t(jj))
          enddo
        endif
      enddo

      if(idb.gt.0) then

        print '(a/(10es11.3))', 'rh_out  ',
     &           (wh(i)/satden(tzero/t(i)),i=1,nz)
        print '(a/(10es11.3))', 'wh_out  ',
     &           (wh(i),i=1,nz)
      endif

      return

      contains

      real(kr) function satden(a)
      real(kr) a
      satden=a*exp(18.916758_kr-a*(14.845878_kr+a*2.4918766_kr))
      return
      end function

      end

c=======================================================================
      subroutine saturate(nz,ncldz,lcld,z,t,rhcld,wh,idb)
      use params, only: kr, tzero
      implicit none
      integer :: lbot, ltop, ncldz, j, idb, i, nz, jj
      real(kr) :: t(*),z(*),wh(*),lcld(*),
     &     cldfac, clrfac, den1, den2, du, dz, tst1, tst2, wvpn,
     &     a, wvp, zbot, ztop, rhcld, wvpclr, wvpcld

c
c purpose: modify the watervapor density inside clouds to have a relative
c          humidity of RHCLD.  reduce water vapor density outside
c          of cloud to maintain total water vapor path
c
c          NOTE: due to normalization procedure the final water vapor
c                density inside the cloud may not exactly match 
c                (rhcld)*(saturated water vapor density), 
c                but it should be close.
c
c INPUT:
c   nz        number of atmospheric levels
c   ncldz     size cloud layer arrays
c   lcld      cloud layer array
c   z         altitude (km) of atmospheric layers (z(nz)=surface)
c   t         temperature at altitude z
c   rhcld     relative humidity of cloud layer, used to specify
c             water vapor density.  
c                                   
c
c h2osat is the mass density (g/m3) of water vapor at 100% saturation 
c (source: handbook of chemistry and physics, h2o vapor pressure table d-112,
c assuming density related to pressure by ideal gas law)
c
c
c            satden(a)=a*exp(18.916758-a*(14.845878+a*2.4918766))
c
c  adjust cloud layer, don't adjust clear layers
c
c  adjust cloud layer and adjust clear layer to compensate for increased wvp
c
      if(idb.gt.0) print '(a/(10es11.3))', 'wh_in  ',
     &           (wh(i)/satden(tzero/t(i)),i=1,nz)

      wvp=0.
      zbot=z(1)
      do i=1,nz
        if(i.eq.1) then
          ztop=.5*(z(2)+zbot)
        elseif(i.eq.nz) then
          ztop=z(nz)
        else
          ztop=.5*(z(i+1)+z(i))
        endif
        wvp=wvp+.1*(ztop-zbot)*wh(i)
        zbot=ztop
      enddo

      if(wvp.eq.0.) then
        print *,'Error in saturate --- ' //
     &       ' original column water vapor is zero -- can not modify'
        stop
      endif

      do i=1,ncldz
        call levrng(ncldz,lcld,i,lbot,ltop)
        if(lbot.ne.0) then
          ! saturate to top of cloud (ltop-1) 
          ! tag saturated points with negative water vapor density
          do j=max(ltop-1,1),lbot ! saturate to top of cloud
            jj=nz-j+1
            wh(jj)=-rhcld*satden(tzero/t(jj)) 
          enddo
        endif
      enddo

      wvpclr=0.
      wvpcld=0.
      zbot=z(1)
      do i=1,nz-1
        if(i.eq.1) then
          ztop=.5*(z(2)+zbot)
        elseif(i.eq.nz) then
          ztop=z(nz)
        else
          ztop=.5*(z(i+1)+z(i))
        endif
        if(wh(i).gt.0.) wvpclr=wvpclr+.1*(ztop-zbot)*wh(i)
        if(wh(i).lt.0.) wvpcld=wvpcld-.1*(ztop-zbot)*wh(i)
        zbot=ztop
      enddo

      if(wvpcld.eq.0) then
        print *,'Error in saturate --- ' //
     &       'water vapor density in cloud = 0 ?'
        stop
      endif

      if(wvpclr.eq.0) then
        cldfac=wvp/wvpcld
        clrfac=1.e-30
      else
        clrfac=(wvp-wvpcld)/wvpclr
        cldfac=1.
        if(clrfac.lt.0) then
          clrfac=1.e-30
          cldfac=wvp/wvpcld
        endif
      endif

      do i=1,nz
        if(wh(i).lt.0) then
          wh(i)=-cldfac*wh(i)
        else
          wh(i)=clrfac*wh(i)
        endif
      enddo

      wvpn=0.
      do i=1,nz-1
        dz=(z(i+1)-z(i))
        den1=wh(i)
        den2=wh(i+1)
        tst1=abs(den1-den2)
        tst2=min(den1,den2)
        if(tst1.le..001*den1 .or. tst2.eq.0.) then
          du=.5*dz*(den1+den2)
        else
          du=dz*(den1-den2)/log(den1/den2)
        endif
        wvpn=wvpn+.1*du
      enddo

      wh(1:nz)=wh(1:nz)*wvp/wvpn

      if(idb.gt.0) then
      
        print '(4a10/4f10.3)', 'wvp','wvpn','wvpclr','wvpcld',
     &       wvp,wvpn,wvpclr,wvpcld
        print '(a/(10es11.3))', 'rh_out  ',
     &           (wh(i)/satden(tzero/t(i)),i=1,nz)
        print '(a/(10es11.3))', 'wh_out  ',
     &           (wh(i),i=1,nz)
      endif      
      return

      contains
      real(kr) function satden(a)
      real(kr) a
      satden=a*exp(18.916758_kr-a*(14.845878_kr+a*2.4918766_kr))
      return
      end function

      end
c===========================================================================
      subroutine modatm(nz,sclh2o,uw,uo3,o3trp,ztrp,pbar,z,p,wh,wo)
c
c purpose:  modify the atmospheric profiles of water vapor, ozone and
c           pressure.  the vertical water vapor and ozone profiles
c           will be scaled by a constant factor so that the integrated
c           amounts will match those specified in the input.  
c
c           when a water vapor scale height is specified the original
c           model atmosphere vertical profile is ignored and the water
c           vapor is distributed with the specified scale hieght and
c           total water vapor content.
c
c           the pressure profile is scaled by the ratio of the
c           specified surface pressure and default surface pressure,
c           p(1), for the model atmosphere. (if the original
c           atmosphere is in pressure equilibrium then the scaled one
c           will be also.)
c
c           Certain input values for SCLH2O, UW, UO3 and PBAR can be used
c           to select default atmospheric profiles:
c 
c           sclh2o <  0     use default water vapor scale height
c           uw     <  0     use default water vapor content
c           uo3    <  0     use default ozone content for z.ge.ztrp
c           o3trp  <  0     use default ozone content for z.lt.ztrp
c           ztrp            altitude of tropopause for ozone adjustment 
c                           (since the default value of ztrp is zero,
c                            uo3 usually controls the total ozone amount)
c           pbar   <= 0     use default surface pressure
c
c
c input:    units
c
c   sclh2o    km        water vapor scale height 
c
c   uw       g/cm**2    water vapor column depth        
c
c   uo3      atm-cm     ozone column depth for z.ge.ztrp. 
c                       ozone density is adjusted up or down to match uo3.  
c
c   o3trp    atm-cm     ozone column depth for z.lt.ztrp.
c
c   ztrp     km         altitude of tropopause for ozone density adjustment
c                       if o3trp.eq.0, uo3 control total ozone amount.
c
c   pbar     mb         surface pressure 
c
c   nz                  number of atmospheric levels
c
c           note: 1 atm-cm = 1 loschmidt / cm**2 = 1000 dobson units
c           
      use params, only: kr, pmo, alosch
      implicit none

      integer :: i, nz

      real(kr) :: z(*),p(*),wh(*),wo(*),
     &   den1, den2, du, dz, fach2o, facp, facstr, factrp, o3trp, 
     &   ofac, pbar, sclh2o, strto3, toth2o, tropo3, tropo3p, uo3, 
     &   uw, w0, ztrp, tst1, tst2

      logical :: debug=.false.
c
c modify water vapor
c
      if(uw.ge.0.) then
c      if(uw.ne.0.) then
        if(debug) write(12,'(a/(10es11.3))') 'h2o before',wh(1:nz)
        if(sclh2o.gt.0.) then
          w0=uw/sclh2o
          do 10 i=1,nz
            wh(i)=w0*exp(-z(i)/sclh2o)
 10       continue
        else
          toth2o=0.
          do 20 i=nz-1,1,-1
            dz=z(i+1)-z(i)
            den1=wh(i)
            den2=wh(i+1)
            tst1=abs(den1-den2)
            tst2=min(den1,den2)
            if(tst1.le..001*den1 .or. tst2.eq.0.) then
              du=.5*dz*(den1+den2)
            else
              du=dz*(den1-den2)/log(den1/den2)
            endif
            toth2o=toth2o+du
 20       continue
c         convert h2o sum from g-km/m3 to g/cm2 
          toth2o=0.1*toth2o
          fach2o=uw/toth2o
          do 30 i=1,nz
            wh(i)=fach2o*wh(i)
 30       continue
        endif
        if(debug) write(12,'(a/(10es11.3))') 'h2o after',wh(1:nz)
      endif
c
c modify ozone
c
      if(uo3.ge.0..or.o3trp.ge.0.) then
        ofac=.5/(30*pmo*alosch)
        if(debug) then    
          write(*,'(a/(10es11.3))') 'z',z(1:nz)
          write(*,'(a/(10es11.3))') 'ozone before',wo(1:nz)
        endif
        tropo3=0.
        tropo3p=0.
        strto3=0.
        do 40 i=nz-1,1,-1
          if(z(i).ge.ztrp) then
            strto3=strto3+ofac*(z(i+1)-z(i))*(wo(i)+wo(i+1))
          elseif(tropo3p.eq.0) then
            tropo3p=ofac*(z(i+1)-z(i))*wo(i+1)
            tropo3=ofac*(z(i+1)-z(i))*wo(i)
          else
            tropo3=tropo3+ofac*(z(i+1)-z(i))*(wo(i)+wo(i+1))
          endif
 40     continue
c         convert o3 sum from g-km/m3 to atm-cm ( = loschmidts-cm)

        if(debug) then
          write(*,'(4a11)') 'tropo3','strto3','total','tropo3p'
          write(*,'(4es11.3)')  tropo3+tropo3p,strto3,
     &       (strto3+tropo3+tropo3p),tropo3p
        endif

        factrp=1.
        facstr=1.
        if(uo3.ge.0.) then
          if(strto3.eq.0) then
            print *, 'Error in modatm -- ' //
     &           'original ozone column density above ZTRP = 0 ' //
     &           ' -- can not modify'
            stop
          endif
          facstr=uo3/strto3
        endif
        if(o3trp.ge.0) then
          if(tropo3.eq.0) then 
            print *,'Error in modatm -- ' //
     &           'original ozone column density below ZTRP = 0 ' //
     &           ' -- can not modify'
            stop
          endif

          ! tropo3*factrp + tropo3p = o3trp
          ! tropo3p represents the column density contribution due to
          ! the wo value just at or above ztrp (i.e., where the factrp
          ! factor has no effect)

          factrp=max(o3trp-tropo3p*facstr,0._kr)/tropo3
          if(debug.and.o3trp.lt.tropo3p*facstr) print '(a)',
     &         'Warning --  (o3trp < tropo3p),  ' //
     &         'full tropospheric ozone adjustment not possible'

        endif

        do i=1,nz
          if(z(i).lt.ztrp) then
            wo(i)=factrp*wo(i)
          else
            wo(i)=facstr*wo(i)
          endif
        enddo

        strto3=0.
        tropo3=0.
        do 50 i=1,nz-1
          if(z(i).ge.ztrp) then
            strto3=strto3+ofac*(z(i+1)-z(i))*(wo(i)+wo(i+1))
          else
            tropo3=tropo3+ofac*(z(i+1)-z(i))*(wo(i)+wo(i+1))
          endif
 50     continue
        if(debug) then
          write(*,'(a/(10es11.3))') 'ozone after',wo(1:nz)
          write(*,'(4a11)') 'tropo3','strto3','total','tropo3p'
          write(*,'(4es11.3)')  tropo3,strto3,
     &       (strto3+tropo3),facstr*tropo3p
          write(*,'(3es11.3,a)') o3trp,uo3,o3trp+uo3,' input values'
        endif
      endif
c
c modify pressure
c
      if(pbar.ge.0.) then
        if(debug) write(*,'(a/(10es11.3))') 'pressure before:',p(1:nz)
        facp=pbar/p(1)
        do 70 i=1,nz
          p(i)=facp*p(i)
 70     continue
        if(debug) write(*,'(a/(10es11.3))') 'pressure after:',p(1:nz)
      endif
      if(debug) stop 'stop in modatm'

      return
      end
C=======================================================================
      subroutine atms(idatm,amix,nz,z,p,t,wh,wo)
c
c     input: 
c       idatm        atmospheric model index, negative to print and quit
c       amix         mixing fraction of user specified atmosphere
c       nz           number of layers
c     output:
c       z            layer altitude
c       p            pressure
c       t            temperature
c       wh           water vapor
c       wo           ozone
c
      use params, only: mxly, kr

      implicit none
      integer :: i, ia, idatm, ierr, nz
      real(kr) :: amix,zz(mxly),pp(mxly),tt(mxly),hh(mxly),oo(mxly),
     &     z(*),p(*),t(*),wh(*),wo(*)

      ia=abs(idatm)

      if(ia.eq.0) call useratm(nz,z,p,t,wh,wo)
      if(ia.eq.1) call tropic(nz,z,p,t,wh,wo)
      if(ia.eq.2) call midsum(nz,z,p,t,wh,wo)
      if(ia.eq.3) call midwin(nz,z,p,t,wh,wo)
      if(ia.eq.4) call subsum(nz,z,p,t,wh,wo)
      if(ia.eq.5) call subwin(nz,z,p,t,wh,wo)
      if(ia.eq.6) call us62(nz,z,p,t,wh,wo)

      if(amix.gt.-1.) then
        call useratm(nz,zz,pp,tt,hh,oo)
        ierr=0
        do 10 i=1,nz
          if(abs(zz(i)-z(i)).gt.0.01) ierr=1
          p(i)=p(i)*(1.-amix)+pp(i)*amix
          t(i)=t(i)*(1.-amix)+tt(i)*amix
          wh(i)=wh(i)*(1.-amix)+hh(i)*amix
          wo(i)=wo(i)*(1.-amix)+oo(i)*amix
 10     continue
        if(ierr.eq.1) stop 'atms -- vertical grids do not match'
      endif

      return
      end
c=======================================================================
      subroutine useratm(nz,z,p,t,wh,wo)

      use params, only: mxly, kr

      implicit none
      integer :: i, nz
      real(kr) :: z(*),p(*),t(*),wh(*),wo(*)             

      open(unit=13,status='old',form='formatted',file='atms.dat')

      read(13,*) nz

      if(nz.gt.mxly) then
        write(*,'(a,i3,a,i3)') 'error in USERATM: ',nz,
     &     ' layers specified in ATMS.DAT, but current limit is ',mxly
        stop
      endif

      do 10 i=nz,1,-1
        read(13,*) z(i),p(i),t(i),wh(i),wo(i)
 10   continue
      close(unit=13)

c     atmosphere should be entered from top to bottom, if not then reverse

      if(z(1).gt.z(nz)) then 
        call reverse(nz,z)
        call reverse(nz,p)
        call reverse(nz,t)
        call reverse(nz,wh)
        call reverse(nz,wo)
      endif

      return
      end
c=========================================================================

      subroutine zgrid(nz,z,p,t,wh,wo,zgrid1,zgrid2,ngrid)

c     module:    zgrid
c
c     purpose:   change vertical resolution of model atmosphere
c
c     input:
c       zgrid1   variable resolution parameter, a floating point number
c                that specifies the grid resolution of the lower part
c                of the grid in kilometers. 
c
c       zgrid2   variable resolution parameter, a floating point number
c                that specifies the maximum allowable grid spacing.
c                For example ZGRID2=40 indicates that the top layer
c                should be 40 km thick and that all other layers should
c                be smaller. 
c
c                A smooth function is used to set the vertical
c                distribution of grid points.  Since this distribution
c                varies smoothly with height, problems associated sharp
c                transitions in step size are avoided (such as
c                spuriously large flux divergence at the jump point).
c
c                if either ZGRID1 or ZGRID2 are specified as a negative
c                numbers REGRID will use the absolute value of the
c                parameters to compute a vertical scale, print the
c                regrided model atmosphere and stop execution.  This
c                option can be used to preview the effect of a given
c                values of ZGRID1 and ZGRID2. 
c
c       ngrid    desired number of vertical grid points.
c
c
c     input/output:
c
c       nz       input as number of vertical grid points in orginal
c                atmosphere, set to NGRID on output 
c
c       z        input:  original vertical position of grid points
c                output: regrided vertical positions
c
c       p        pressure at height z
c       t        temperature at heigth z
c       wh       water vapor density at height z
c       wo       ozone density at height z
c
c       
c
      use params, only: mxly, zero, one, kr

      implicit none
      integer :: i, j, jj, ngrid, nz, ng
      real(kr) :: z(*),p(*),t(*),wh(*),wo(*), 
     &     a, beta, fz, toprat, x, zgrid1, zgrid2, ztop,
     &     zz(mxly),pp(mxly),tt(mxly),whh(mxly),woo(mxly)
     
      real(kr), parameter :: tol=0.99

      ng=min(mxly,abs(ngrid))

      a=zgrid1*float(ng-1)
      ztop=z(nz)
      if( a .ge. ztop) then 
        a=ztop
        beta=0.
      else
        a=min(a,tol*(ztop-zgrid2))
        toprat=float(ng-2)/(ng-1)
        beta=log(((ztop-zgrid2)/(toprat*a)-1.)/
     &       (ztop/a-1.))/log(toprat)
      endif

      j=2

      do 30 i=1,ng
        x=float(i-1)/float(ng-1)
        zz(i)=a*x*(1.+(ztop/a-1.)*x**beta)
        do 10 jj=j,nz
          if(zz(i).le.z(jj)) goto 20
 10     continue
 20     continue
        j=min(jj,nz)
        fz=(zz(i)-z(j-1))/(z(j)-z(j-1))
        fz=min(max(fz,zero),one)
        
c         write(*,'(2i10,3f12.4)') i,j,zz(i),z(j),fz
        
        pp(i)  =  p(j-1)  * (p(j)/p(j-1))**fz
        tt(i)  =  t(j-1)  * (1.-fz) + t(j)  * fz
        if(min(wh(j),wh(j-1)).gt.0.) then
          whh(i)= wh(j-1) * (wh(j)/wh(j-1))** fz
        else
          whh(i) =  wh(j-1) * (1.-fz) + wh(j) * fz
        endif
        if(min(wo(j),wo(j-1)).gt.0.) then
          woo(i)= wo(j-1) * (wo(j)/wo(j-1))** fz
        else
          woo(i) =  wo(j-1) * (1.-fz) + wo(j) * fz
        endif
 30   continue
      
      do 40 i=1,ng
        z(i)=zz(i)
        p(i)=pp(i)
        t(i)=tt(i)
        wh(i)=whh(i)
        wo(i)=woo(i)
 40   continue
      
      nz=ng
      
      return
      end

c=========================================================================
      subroutine elevate(nz,z,p,t,wh,wo,zsurf)

c     module:    elevate
c
c     purpose:   change surface altitude of model atmosphere if the new
c                altitude of the surface is between the two points
c                original in the original atmospheric grid, values for
c                that altitude are obtained by interpolation.  Thus the
c                grid spacing of of the lowest layer may be made
c                smaller. All points above the surface grid are
c                renumbered.
c                
c
c     input:
c       zsurf    altitude of the surface layer.  
c
c     input/output:
c       nz       input as number of vertical grid points in orginal
c                atmosphere, set to NGRID on output 
c
c       z        input:  original vertical position of grid points
c                output: regrided or renumbered vertical positions
c
c       p        pressure at height z
c       t        temperature at heigth z
c       wh       water vapor density at height z
c       wo       ozone density at height z
c
c       
c
      use params, only: mxly, kr

      implicit none
      real(kr), parameter :: eps=.001

      integer :: i, j, ks, kz, nz

      real(kr) :: dist, dist1, dist2, fz, zsurf,
     &     z(*),p(*),t(*),wh(*),wo(*),
     &     zz(mxly),pp(mxly),tt(mxly),whh(mxly),woo(mxly)
     
      call locate(z,nz,zsurf,kz)
      dist1=zsurf-z(kz)
      dist2=z(kz+1)-zsurf
      dist=z(kz+1)-z(kz)
      if(dist1.lt.eps*dist) then
        zz(1)=z(kz)
        fz=0.
        ks=kz+1
      elseif(dist2.lt.eps*dist) then 
        zz(1)=z(kz+1)
        fz=1.
        ks=kz+2
      else
        zz(1)=zsurf
        fz=(zsurf-z(kz))/(z(kz+1)-z(kz))
        ks=kz+1
      endif
      pp(1)  =  p(kz)  * (p(kz+1)/p(kz))**fz
      tt(1)  =  t(kz)  * (1.-fz) + t(kz+1)  * fz

      if(min(wh(kz+1),wh(kz)).gt.0.) then
        whh(1)= wh(kz) * (wh(kz+1)/wh(kz))** fz
      else
        whh(1) =  wh(kz) * (1.-fz) + wh(kz+1) * fz
      endif
      if(min(wo(kz+1),wo(kz)).gt.0.) then
        woo(1)= wo(kz) * (wo(kz+1)/wo(kz))** fz
      else
        woo(1) =  wo(kz) * (1.-fz) + wo(kz+1) * fz
      endif


      i=1
      do 30 j=ks,nz
        i=i+1
        zz(i)=z(j)
        pp(i)=p(j)
        tt(i)=t(j)
        whh(i)=wh(j)
        woo(i)=wo(j)
 30   continue

      nz=i

      do 40 i=1,nz
        z(i)=zz(i)
        p(i)=pp(i)
        t(i)=tt(i)
        wh(i)=whh(i)
        wo(i)=woo(i)
 40   continue


      return
      end
c=========================================================================
      subroutine tropic(nz,z,p,t,wh,wo)
      use params, only: kr
      implicit none
      integer, parameter :: mz=33

      integer :: i, nz

      real(kr) :: z(*),p(*),t(*),wh(*),wo(*),
     &     z1(mz),p1(mz),t1(mz),wh1(mz),wo1(mz)
c                                                                       
c     model: tropical mc'clatchey                                       
c                                                                       
      data z1/
     &    0.,    1.,    2.,    3.,    4.,    5.,    6.,    7.,    8.,   
     &    9.,   10.,   11.,   12.,   13.,   14.,   15.,   16.,   17.,   
     &   18.,   19.,   20.,   21.,   22.,   23.,   24.,   25.,   30.,   
     &   35.,   40.,   45.,   50.,   70.,  100./                 
      data p1/
     &  1.013e+03,9.040e+02,8.050e+02,7.150e+02,6.330e+02,5.590e+02,      
     &  4.920e+02,4.320e+02,3.780e+02,3.290e+02,2.860e+02,2.470e+02,      
     &  2.130e+02,1.820e+02,1.560e+02,1.320e+02,1.110e+02,9.370e+01,      
     &  7.890e+01,6.660e+01,5.650e+01,4.800e+01,4.090e+01,3.500e+01,      
     &  3.000e+01,2.570e+01,1.220e+01,6.000e+00,3.050e+00,1.590e+00,      
     &  8.540e-01,5.790e-02,3.000e-04/                          
      data t1/
     &  3.000e+02,2.940e+02,2.880e+02,2.840e+02,2.770e+02,2.700e+02,      
     &  2.640e+02,2.570e+02,2.500e+02,2.440e+02,2.370e+02,2.300e+02,      
     &  2.240e+02,2.170e+02,2.100e+02,2.040e+02,1.970e+02,1.950e+02,      
     &  1.990e+02,2.030e+02,2.070e+02,2.110e+02,2.150e+02,2.170e+02,      
     &  2.190e+02,2.210e+02,2.320e+02,2.430e+02,2.540e+02,2.650e+02,      
     &  2.700e+02,2.190e+02,2.100e+02/                          
      data wh1/
     &  1.900e+01,1.300e+01,9.300e+00,4.700e+00,2.200e+00,1.500e+00,      
     &  8.500e-01,4.700e-01,2.500e-01,1.200e-01,5.000e-02,1.700e-02,      
     &  6.000e-03,1.800e-03,1.000e-03,7.600e-04,6.400e-04,5.600e-04,      
     &  5.000e-04,4.900e-04,4.500e-04,5.100e-04,5.100e-04,5.400e-04,      
     &  6.000e-04,6.700e-04,3.600e-04,1.100e-04,4.300e-05,1.900e-05,      
     &  6.300e-06,1.400e-07,1.000e-09/                          
      data wo1/
     &  5.600e-05,5.600e-05,5.400e-05,5.100e-05,4.700e-05,4.500e-05,      
     &  4.300e-05,4.100e-05,3.900e-05,3.900e-05,3.900e-05,4.100e-05,      
     &  4.300e-05,4.500e-05,4.500e-05,4.700e-05,4.700e-05,6.900e-05,      
     &  9.000e-05,1.400e-04,1.900e-04,2.400e-04,2.800e-04,3.200e-04,      
     &  3.400e-04,3.400e-04,2.400e-04,9.200e-05,4.100e-05,1.300e-05,      
     &  4.300e-06,8.600e-08,4.300e-11/                          
      do 1 i=1,mz                                                       
        z(i)=z1(i)                                                        
        p(i)=p1(i)                                                        
        t(i)=t1(i)                                                        
        wh(i)=wh1(i)                                                      
        wo(i)=wo1(i)                                                      
    1 continue  
      nz=mz                                                        
      return                                                            
      end                                                               
c=======================================================================
      subroutine midsum(nz,z,p,t,wh,wo)
      use params, only: kr
      implicit none
      integer, parameter :: mz=33
      
      integer :: i, nz

      real(kr) :: z(*),p(*),t(*),wh(*),wo(*),
     &     z2(mz),p2(mz),t2(mz),wh2(mz),wo2(mz)
c                                                                       
c     model: midlatitude summer mc'clatchey                             
c                                                                       
      data z2/
     &    0.,    1.,    2.,    3.,    4.,    5.,    6.,    7.,    8.,   
     &    9.,   10.,   11.,   12.,   13.,   14.,   15.,   16.,   17.,   
     &   18.,   19.,   20.,   21.,   22.,   23.,   24.,   25.,   30.,   
     &   35.,   40.,   45.,   50.,   70.,  100./                 
      data p2/
     &  1.013e+03,9.020e+02,8.020e+02,7.100e+02,6.280e+02,5.540e+02,      
     &  4.870e+02,4.260e+02,3.720e+02,3.240e+02,2.810e+02,2.430e+02,      
     &  2.090e+02,1.790e+02,1.530e+02,1.300e+02,1.110e+02,9.500e+01,      
     &  8.120e+01,6.950e+01,5.950e+01,5.100e+01,4.370e+01,3.760e+01,      
     &  3.220e+01,2.770e+01,1.320e+01,6.520e+00,3.330e+00,1.760e+00,      
     &  9.510e-01,6.710e-02,3.000e-04/                          
      data t2/
     &  2.940e+02,2.900e+02,2.850e+02,2.790e+02,2.730e+02,2.670e+02,      
     &  2.610e+02,2.550e+02,2.480e+02,2.420e+02,2.350e+02,2.290e+02,      
     &  2.220e+02,2.160e+02,2.160e+02,2.160e+02,2.160e+02,2.160e+02,      
     &  2.160e+02,2.170e+02,2.180e+02,2.190e+02,2.200e+02,2.220e+02,      
     &  2.230e+02,2.240e+02,2.340e+02,2.450e+02,2.580e+02,2.700e+02,      
     &  2.760e+02,2.180e+02,2.100e+02/                          
      data wh2/
     &  1.400e+01,9.300e+00,5.900e+00,3.300e+00,1.900e+00,1.000e+00,      
     &  6.100e-01,3.700e-01,2.100e-01,1.200e-01,6.400e-02,2.200e-02,      
     &  6.000e-03,1.800e-03,1.000e-03,7.600e-04,6.400e-04,5.600e-04,      
     &  5.000e-04,4.900e-04,4.500e-04,5.100e-04,5.100e-04,5.400e-04,      
     &  6.000e-04,6.700e-04,3.600e-04,1.100e-04,4.300e-05,1.900e-05,      
     &  1.300e-06,1.400e-07,1.000e-09/                          
      data wo2/
     &  6.000e-05,6.000e-05,6.000e-05,6.200e-05,6.400e-05,6.600e-05,      
     &  6.900e-05,7.500e-05,7.900e-05,8.600e-05,9.000e-05,1.100e-04,      
     &  1.200e-04,1.500e-04,1.800e-04,1.900e-04,2.100e-04,2.400e-04,      
     &  2.800e-04,3.200e-04,3.400e-04,3.600e-04,3.600e-04,3.400e-04,      
     &  3.200e-04,3.000e-04,2.000e-04,9.200e-05,4.100e-05,1.300e-05,      
     &  4.300e-06,8.600e-08,4.300e-11/                          

      do 1 i=1,mz                                                       
        z(i)=z2(i)                                                        
        p(i)=p2(i)                                                        
        t(i)=t2(i)                                                        
        wh(i)=wh2(i)                                                      
        wo(i)=wo2(i)                                                      
    1 continue             
      nz=mz                                             
      return                                                            
      end                                                               
c=========================================================================
      subroutine midwin(nz,z,p,t,wh,wo)
      use params, only: kr
      implicit none
      integer, parameter :: mz=33
      
      integer :: i, nz

      real(kr) :: z(*),p(*),t(*),wh(*),wo(*),
     &     z3(mz),p3(mz),t3(mz),wh3(mz),wo3(mz)
c                                                                       
c     model: midlatitude winter mc'clatchey                             
c                                                                       
      data z3/
     &     0.,    1.,    2.,    3.,    4.,    5.,    6.,    7.,    8.,   
     &     9.,   10.,   11.,   12.,   13.,   14.,   15.,   16.,   17.,   
     &    18.,   19.,   20.,   21.,   22.,   23.,   24.,   25.,   30.,   
     &    35.,   40.,   45.,   50.,   70.,  100./                 
      data p3/
     &     1.018e+03,8.973e+02,7.897e+02,6.938e+02,6.081e+02,5.313e+02,      
     &     4.627e+02,4.016e+02,3.473e+02,2.992e+02,2.568e+02,2.199e+02,      
     &     1.882e+02,1.610e+02,1.378e+02,1.178e+02,1.007e+02,8.610e+01,      
     &     7.350e+01,6.280e+01,5.370e+01,4.580e+01,3.910e+01,3.340e+01,      
     &     2.860e+01,2.430e+01,1.110e+01,5.180e+00,2.530e+00,1.290e+00,      
     &  6.820e-01,4.670e-02,3.000e-04/                          
      data t3/
     &  2.722e+02,2.687e+02,2.652e+02,2.617e+02,2.557e+02,2.497e+02,      
     &  2.437e+02,2.377e+02,2.317e+02,2.257e+02,2.197e+02,2.192e+02,      
     &  2.187e+02,2.182e+02,2.177e+02,2.172e+02,2.167e+02,2.162e+02,      
     &  2.157e+02,2.152e+02,2.152e+02,2.152e+02,2.152e+02,2.152e+02,      
     &  2.152e+02,2.152e+02,2.174e+02,2.278e+02,2.432e+02,2.585e+02,      
     &  2.657e+02,2.307e+02,2.102e+02/                          
      data wh3/
     &  3.500e+00,2.500e+00,1.800e+00,1.200e+00,6.600e-01,3.800e-01,      
     &  2.100e-01,8.500e-02,3.500e-02,1.600e-02,7.500e-03,6.900e-03,      
     &  6.000e-03,1.800e-03,1.000e-03,7.600e-04,6.400e-04,5.600e-04,      
     &  5.000e-04,4.900e-04,4.500e-04,5.100e-04,5.100e-04,5.400e-04,      
     &  6.000e-04,6.700e-04,3.600e-04,1.100e-04,4.300e-05,1.900e-05,      
     &  6.300e-06,1.400e-07,1.000e-09/                          
      data wo3/
     &  6.000e-05,5.400e-05,4.900e-05,4.900e-05,4.900e-05,5.800e-05,      
     &  6.400e-05,7.700e-05,9.000e-05,1.200e-04,1.600e-04,2.100e-04,      
     &  2.600e-04,3.000e-04,3.200e-04,3.400e-04,3.600e-04,3.900e-04,      
     &  4.100e-04,4.300e-04,4.500e-04,4.300e-04,4.300e-04,3.900e-04,      
     &  3.600e-04,3.400e-04,1.900e-04,9.200e-05,4.100e-05,1.300e-05,      
     &  4.300e-06,8.600e-08,4.300e-11/                          
      do 1 i=1,mz                                                       
        z(i)=z3(i)                                                        
        p(i)=p3(i)                                                        
        t(i)=t3(i)                                                        
        wh(i)=wh3(i)                                                      
        wo(i)=wo3(i)                                                      
    1 continue             
      nz=mz                                             
      return                                                            
      end                                                               
c=========================================================================
      subroutine subsum(nz,z,p,t,wh,wo)
      use params, only: kr
      implicit none
      integer, parameter :: mz=33
      
      integer :: i, nz

      real(kr) :: z(*),p(*),t(*),wh(*),wo(*),
     &     z4(mz),p4(mz),t4(mz),wh4(mz),wo4(mz)
c                                                                       
c     model: subarctic summer mc'clatchey                             
c                                                                       
      data z4/
     &    0.,    1.,    2.,    3.,    4.,    5.,    6.,    7.,    8.,   
     &    9.,   10.,   11.,   12.,   13.,   14.,   15.,   16.,   17.,   
     &   18.,   19.,   20.,   21.,   22.,   23.,   24.,   25.,   30.,   
     &   35.,   40.,   45.,   50.,   70.,  100./                 
      data p4/
     &  1.010e+03,8.960e+02,7.929e+02,7.000e+02,6.160e+02,5.410e+02,      
     &  4.730e+02,4.130e+02,3.590e+02,3.107e+02,2.677e+02,2.300e+02,      
     &  1.977e+02,1.700e+02,1.460e+02,1.250e+02,1.080e+02,9.280e+01,      
     &  7.980e+01,6.860e+01,5.890e+01,5.070e+01,4.360e+01,3.750e+01,      
     &  3.227e+01,2.780e+01,1.340e+01,6.610e+00,3.400e+00,1.810e+00,      
     &  9.870e-01,7.070e-02,3.000e-04/                          
      data t4/
     &  2.870e+02,2.820e+02,2.760e+02,2.710e+02,2.660e+02,2.600e+02,      
     &  2.530e+02,2.460e+02,2.390e+02,2.320e+02,2.250e+02,2.250e+02,      
     &  2.250e+02,2.250e+02,2.250e+02,2.250e+02,2.250e+02,2.250e+02,      
     &  2.250e+02,2.250e+02,2.250e+02,2.250e+02,2.250e+02,2.250e+02,      
     &  2.260e+02,2.280e+02,2.350e+02,2.470e+02,2.620e+02,2.740e+02,      
     &  2.770e+02,2.160e+02,2.100e+02/                          
      data wh4/
     &  9.100e+00,6.000e+00,4.200e+00,2.700e+00,1.700e+00,1.000e+00,      
     &  5.400e-01,2.900e-01,1.300e-01,4.200e-02,1.500e-02,9.400e-03,      
     &  6.000e-03,1.800e-03,1.000e-03,7.600e-04,6.400e-04,5.600e-04,      
     &  5.000e-04,4.900e-04,4.500e-04,5.100e-04,5.100e-04,5.400e-04,      
     &  6.000e-04,6.700e-04,3.600e-04,1.100e-04,4.300e-05,1.900e-05,      
     &  6.300e-06,1.400e-07,1.000e-09/                          
      data wo4/
     &  4.900e-05,5.400e-05,5.600e-05,5.800e-05,6.000e-05,6.400e-05,      
     &  7.100e-05,7.500e-05,7.900e-05,1.100e-04,1.300e-04,1.800e-04,      
     &  2.100e-04,2.600e-04,2.800e-04,3.200e-04,3.400e-04,3.900e-04,      
     &  4.100e-04,4.100e-04,3.900e-04,3.600e-04,3.200e-04,3.000e-04,      
     &  2.800e-04,2.600e-04,1.400e-04,9.200e-05,4.100e-05,1.300e-05,      
     &  4.300e-06,8.600e-08,4.300e-11/                          
      do 1 i=1,mz                                                       
        z(i)=z4(i)                                                        
        p(i)=p4(i)                                                        
        t(i)=t4(i)                                                        
        wh(i)=wh4(i)                                                      
        wo(i)=wo4(i)                                                      
    1 continue             
      nz=mz                                             
      return                                                            
      end                                                               
c=========================================================================
      subroutine subwin(nz,z,p,t,wh,wo)
      use params, only: kr
      implicit none
      integer, parameter :: mz=33
      
      integer :: i, nz

      real(kr) :: z(*),p(*),t(*),wh(*),wo(*),
     &     z5(mz),p5(mz),t5(mz),wh5(mz),wo5(mz)
c                                                                       
c     model: subarctique winter mc clatchey                             
c                                                                       
      data z5/
     &    0.,    1.,    2.,    3.,    4.,    5.,    6.,    7.,    8.,   
     &    9.,   10.,   11.,   12.,   13.,   14.,   15.,   16.,   17.,   
     &   18.,   19.,   20.,   21.,   22.,   23.,   24.,   25.,   30.,   
     &   35.,   40.,   45.,   50.,   70.,  100./                 
      data p5/
     &  1.013e+03,8.878e+02,7.775e+02,6.798e+02,5.932e+02,5.158e+02,      
     &  4.467e+02,3.853e+02,3.308e+02,2.829e+02,2.418e+02,2.067e+02,      
     &  1.766e+02,1.510e+02,1.291e+02,1.103e+02,9.431e+01,8.058e+01,      
     &  6.882e+01,5.875e+01,5.014e+01,4.277e+01,3.647e+01,3.109e+01,      
     &  2.649e+01,2.256e+01,1.020e+01,4.701e+00,2.243e+00,1.113e+00,      
     &  5.719e-01,4.016e-02,3.000e-04/                          
      data t5/
     &  2.571e+02,2.591e+02,2.559e+02,2.527e+02,2.477e+02,2.409e+02,      
     &  2.341e+02,2.273e+02,2.206e+02,2.172e+02,2.172e+02,2.172e+02,      
     &  2.172e+02,2.172e+02,2.172e+02,2.172e+02,2.166e+02,2.160e+02,      
     &  2.154e+02,2.148e+02,2.141e+02,2.136e+02,2.130e+02,2.124e+02,      
     &  2.118e+02,2.112e+02,2.160e+02,2.222e+02,2.347e+02,2.470e+02,      
     &  2.593e+02,2.457e+02,2.100e+02/                          
      data wh5/
     &  1.200e+00,1.200e+00,9.400e-01,6.800e-01,4.100e-01,2.000e-01,      
     &  9.800e-02,5.400e-02,1.100e-02,8.400e-03,5.500e-03,3.800e-03,      
     &  2.600e-03,1.800e-03,1.000e-03,7.600e-04,6.400e-04,5.600e-04,      
     &  5.000e-04,4.900e-04,4.500e-04,5.100e-04,5.100e-04,5.400e-04,      
     &  6.000e-04,6.700e-04,3.600e-04,1.100e-04,4.300e-05,1.900e-05,      
     &  6.300e-06,1.400e-07,1.000e-09/                          
      data wo5/
     &  4.100e-05,4.100e-05,4.100e-05,4.300e-05,4.500e-05,4.700e-05,      
     &  4.900e-05,7.100e-05,9.000e-05,1.600e-04,2.400e-04,3.200e-04,      
     &  4.300e-04,4.700e-04,4.900e-04,5.600e-04,6.200e-04,6.200e-04,      
     &  6.200e-04,6.000e-04,5.600e-04,5.100e-04,4.700e-04,4.300e-04,      
     &  3.600e-04,3.200e-04,1.500e-04,9.200e-05,4.100e-05,1.300e-05,      
     &  4.300e-06,8.600e-08,4.300e-11/                          
      do 1 i=1,mz                                                       
        z(i)=z5(i)                                                        
        p(i)=p5(i)                                                        
        t(i)=t5(i)                                                        
        wh(i)=wh5(i)                                                      
        wo(i)=wo5(i)                                                      
    1 continue             
      nz=mz                                             
      return                                                            
      end                                                               
c=======================================================================
      subroutine us62(nz,z,p,t,wh,wo)
      use params, only: kr
      implicit none
      integer, parameter :: mz=33
      
      integer :: i, nz

      real(kr) :: z(*),p(*),t(*),wh(*),wo(*),
     &     z6(mz),p6(mz),t6(mz),wh6(mz),wo6(mz)
c                                                                       
c     model: us standard 62 mc clatchey                                 
c                                                                       
      data z6/
     &    0.,    1.,    2.,    3.,    4.,    5.,    6.,    7.,    8.,   
     &    9.,   10.,   11.,   12.,   13.,   14.,   15.,   16.,   17.,   
     &   18.,   19.,   20.,   21.,   22.,   23.,   24.,   25.,   30.,   
     &   35.,   40.,   45.,   50.,   70.,  100./                 
      data p6/
     &  1.013e+03,8.986e+02,7.950e+02,7.012e+02,6.166e+02,5.405e+02,      
     &  4.722e+02,4.111e+02,3.565e+02,3.080e+02,2.650e+02,2.270e+02,      
     &  1.940e+02,1.658e+02,1.417e+02,1.211e+02,1.035e+02,8.850e+01,      
     &  7.565e+01,6.467e+01,5.529e+01,4.729e+01,4.047e+01,3.467e+01,      
     &  2.972e+01,2.549e+01,1.197e+01,5.746e+00,2.871e+00,1.491e+00,      
     &  7.978e-01,5.520e-02,3.008e-04/                          
      data t6/
     &  2.881e+02,2.816e+02,2.751e+02,2.687e+02,2.622e+02,2.557e+02,      
     &  2.492e+02,2.427e+02,2.362e+02,2.297e+02,2.232e+02,2.168e+02,      
     &  2.166e+02,2.166e+02,2.166e+02,2.166e+02,2.166e+02,2.166e+02,      
     &  2.166e+02,2.166e+02,2.166e+02,2.176e+02,2.186e+02,2.196e+02,      
     &  2.206e+02,2.216e+02,2.265e+02,2.365e+02,2.534e+02,2.642e+02,      
     &  2.706e+02,2.197e+02,2.100e+02/                          
      data wh6/
     &  5.900e+00,4.200e+00,2.900e+00,1.800e+00,1.100e+00,6.400e-01,      
     &  3.800e-01,2.100e-01,1.200e-01,4.600e-02,1.800e-02,8.200e-03,      
     &  3.700e-03,1.800e-03,8.400e-04,7.200e-04,6.100e-04,5.200e-04,      
     &  4.400e-04,4.400e-04,4.400e-04,4.800e-04,5.200e-04,5.700e-04,      
     &  6.100e-04,6.600e-04,3.800e-04,1.600e-04,6.700e-05,3.200e-05,      
     &  1.200e-05,1.500e-07,1.000e-09/                          
      data wo6/
     &  5.400e-05,5.400e-05,5.400e-05,5.000e-05,4.600e-05,4.600e-05,      
     &  4.500e-05,4.900e-05,5.200e-05,7.100e-05,9.000e-05,1.300e-04,      
     &  1.600e-04,1.700e-04,1.900e-04,2.100e-04,2.400e-04,2.800e-04,      
     &  3.200e-04,3.500e-04,3.800e-04,3.800e-04,3.900e-04,3.800e-04,      
     &  3.600e-04,3.400e-04,2.000e-04,1.100e-04,4.900e-05,1.700e-05,      
     &  4.000e-06,8.600e-08,4.300e-11/                          
      do 1 i=1,mz                                                       
        z(i)=z6(i)                                                        
        p(i)=p6(i)                                                        
        t(i)=t6(i)                                                        
        wh(i)=wh6(i)                                                      
        wo(i)=wo6(i)                                                      
    1 continue           
      nz=mz                                               
      return                                                            
      end                                                               
c=======================================================================















