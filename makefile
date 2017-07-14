PROG =	sbdart

OBJS =  atms.o disort.o disutil.o drt.o params.o spectra.o \
        tauaero.o taucloud.o taugas.o

FC = /usr/local/bin/gfortran
FFLAGS = -O -C -std=f2003 -Wuninitialized
LDFLAGS =

LIBS =

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

atms.o:  params.o atms.f

disort.o:  params.o disort.f

disutil.o:  params.o disutil.f

drt.o:  params.o tauaero.o taugas.o spectra.o drt.f

params.o:  params.f

spectra.o:  params.o spectra.f

tauaero.o:  params.o tauaero.f

taucloud.o:  params.o taucloud.f

taugas.o:  params.o taugas.f


clean:
	rm -f $(PROG) $(OBJS) *.mod

