      module params
      implicit none

      save   

      integer, parameter :: kr=selected_real_kind(10) 

      integer, parameter ::
     &     mxly=65,             ! maximum number of layers
     &     maxmom=299,          ! maximum number of legendre moments
     &     nstrms=40,           ! maximum number of radiation streams
     &     ncldz=5,             ! number of cloud layers
     &     naerz=5,             ! number of stratospheric aerosol layers
     &     mxq=63,              ! number of lowtran species
     &     mxkd=20              ! max size of k-distribution arrays

      real(kr), parameter ::
     &     pzero=1013.25,       ! standard pressure (mb)
     &     tzero=273.15,        ! standard temperature (K)
     &     re=6371.2,           ! radius of earth (km)
     &     pmo=2.6568e-23,      ! oxygen atomic mass (g)
     &     grav=9.80665,        ! gravitational acceleration (m/s/s)
     &     alosch=2.6868e19,    ! loschmidt number molecules/cm3 at stp
     &     zip=-1.,             ! undefined read value
     &     nan=huge(0.),        ! another undefined read value
     &     one=1.,              ! 1 of selected real kind kr
     &     zero=0.,             ! 0 of selected real kind kr
     &     wl55=0.55,           ! standard visible wavelength
     &     pi=3.1415926536_kr   ! circumference/diameter of a circle 
                                ! the _kr ensures precision matches 
                                ! selected real kind kr

      end module params
