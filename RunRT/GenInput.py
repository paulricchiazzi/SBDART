from collections import OrderedDict
class GenInput:
    def __init__(self):

        self.rtRange = self.setRtRange()
        self.cycles=[] # number of iterations in each for loop
        self.rtvar=[]  # variable rt parameter   [ {"TCLOUD":[0,10,20]}, {"SZA":[0,30,60]}, {"ALBCON":[0,0.5,1.0]}]
        self.rtcons={} # constant rt parameters  { "NRE": 10}
        self.IOUTformat = 10


    def CycleSequence(self):
        """
        given rtvar=[ {"TCLOUD":[0,10,20]}, {"SZA":[0,30,60]}, {"ALBCON":[0,0.5,1.0]}]
        :return: if iout==10: ["SZA","ALBCON"] ,  else: ["TCLOUD", "SZA", "ALBCON"]
        """
        seq=[]
        for hdict in self.rtvar:
            seq.append(hdict.keys()[0])

        return seq[1:] if self.IOUTformat == 10 else seq

    def GetParmMenuItems(self):
        """
        Get the sequence of individual variant values for each nesting cycle directly from command text
        Given a cmd text such as
                         A=0;1
                         B=2;3
                         C=4;5
                         D=6

        :return: [["A=0","A=1"]["B=2","B=3"],["C=4","C=5"],["D=6"]]
        """
        tags=[]
        ncycles = len(self.rtvar)
        if ncycles == 0: return tags
        i=-1
        for cdict in self.rtvar:
            i+=1
            k=cdict.keys()[0]
            tags.append([])
            for v in cdict[k]:
                tags[i].append("{}={}".format(k,v))
        return tags

    def CycleSetup(self, cmd):
        """
        cycle through Runrt command buffer
        :param:cmd: -- RunRT command sequence
        :returns:xvariable: list of values for the first variant if IOUT=10,
        :returns:xlabel: the parameter name used to label x-axis
        """
        lines = cmd.split('\n')
        self.rtvar=[]                   # a list of dictionaries, a dictionary for each nested loop
        self.rtcons={}
        self.cycles=[]
        nesting = -1
        niter = 0
        xvariable = []
        xlabel = ""
        prevcycles = 0                 # check that lock-step variables have same number of iterations
        self.IOUTformat = 10  # this is the default
        for line in lines:
            line = line.split('#')[0].strip()   # strip off comments
            if len(line) == 0:
                continue
            if not line.endswith("&"):
                nesting+=1
            if "=" in line:
                parm, rhs = line.split('=')
                if parm == "IOUT":
                    if rhs.count(';') > 0:
                        rhs=rhs.split(';')[0]                   # only one output format allowed
                    self.IOUTformat = int(rhs)

                if ";" in rhs:
                    covariant = rhs.endswith("&")
                    if covariant:
                        rhs=rhs[0:-2]
                    values=[vv.strip().replace(' ','_') for vv in rhs.split(';')]
                    if covariant and not prevcycles == len(values) :
                        xvariable = values
                        xlabel = "Error: Number of elements in covariant variable, {}, is incorrect.".format(parm)
                        return xvariable, xlabel

                    prevcycles = len(values)
                    #values=[vv.strip() for vv in rhs.split(';')]
                    if not covariant:
                        self.cycles.append(len(values))
                        self.rtvar.append(OrderedDict())
                    self.rtvar[nesting][parm]=values
                    if nesting == 0:
                        # note that xvariable is set here for IOUT=10,
                        # it is reset in Plot01 and Plot11 to wavelength and altitude, respectively
                        xvariable = values[:]
                        xlabel = str(parm)
                else:
                    self.rtcons[parm]=rhs
        return xvariable, xlabel

    def Niter(self):
        """
        :return: total number of iterations required to cycle through command file
        """
        niter = 1
        for n in self.cycles:
            niter*=n
        return niter

    def CycleInput(self, iteration):
        """
        generates INPUT for a given iteration

        :param
            iteration   - iteration number between 0 and Niter()-1
        :returns:
            rtinp        - SBDART input for this iteration
            rtlist       - rtinp broken into list of leading variants in each
                            cycle.  E.g., if cmd is

            TCLOUD=0;10;100
            WLINF=0.5;0.8
            WLSUP=0.5;0.8 &
            ALBCON=0.5;1.0

            rtinp=TCLOUD=0\nWLINF=0.8\nWLSUP=0.8\nALBCON=0.5
            rtlist=["TCLOUD=0","WLINF=0.8","ALBCON=0.5"] for iteration 3
        """
        rtinp=""
        rtlist=[]
        niter = self.Niter()
        if iteration >= niter:
            return None
        else:
            iter = iteration
            for nest in range(0,len(self.rtvar)):
                if nest > 0:
                    iter=iter/self.cycles[nest-1]
                i = iter % self.cycles[nest]
                first = True
                for p,vlist in self.rtvar[nest].iteritems():
                    v=vlist[i]
                    if p[0].isalpha():           # don't write output RT parm unless its legit fortran var
                        rtinp+="{}={}\n".format(p,v)
                    if first:
                        rtlist.append("{}={}".format(p,v))
                    first = False
            for p,v in self.rtcons.iteritems():
                rtinp+="{}={}\n".format(p,v.split('#')[0])

        return rtinp,rtlist

    def DocString(self, key):
        docstr = {"IOUT=1"   : "Spectral output between WLINF and WLSUP",
                  "IOUT=2"   : "Spectral profile of lowtran optical depth",
                  "IOUT=10"  : "Radiant flux at top and bottom of atmosphere",
                  "IOUT=11"  : "Radiant flux at each atmospheric layer",
                  "IOUT=20"  : "Radiance at top of atmosphere",
                  "IOUT=21"  : "Radiance at bottom of atmosphere",
                  "IAER=0"   : "No boundary layer aerosol",
                  "IAER=1"   : "Rural aerosol",
                  "IAER=2"   : "Urban aerosol",
                  "IAER=3"   : "Oceanic aerosol",
                  "IAER=4"   : "Tropospheric aerosol",
                  "IAER=5"   : "User defined (see WLBAER, QBAER, TBAER, WBAER, GBAER)",
                  "JAER=0"   : "No stratospheric aerosol",
                  "JAER=1"   : "Background stratospheric aerosol",
                  "JAER=2"   : "Aged volcanic stratospheric aerosol",
                  "JAER=3"   : "Fresh volcanic stratospheric aerosol",
                  "JAER=4"   : "Meteor dust stratospheric aerosol",
                  "ISALB=-1" : "Read surface albedo from albedo.dat",
                  "ISALB=0"  : "Constant albedo set by ALBCON",
                  "ISALB=1"  : "Snow",
                  "ISALB=2"  : "Clear water",
                  "ISALB=3"  : "Lake water",
                  "ISALB=4"  : "Sea water",
                  "ISALB=5"  : "Sand",
                  "ISALB=6"  : "Vegetation",
                  "ISALB=7"  : "Ocean water BRDF, requires SC parameters",
                  "ISALB=8"  : "Hapke BRDF model, requires SC parameters",
                  "ISALB=9"  : "Ross-thick Li-sparse BRDF, requires SC parameters",
                  "ISALB=10" : "Snow, seawater, sand and vegetation, SC sets partitions",
                  "ISAT=-4"  : "Gaussian filter, WLINF-2*WLSUP to WLINF+2*WLSUP",
                  "ISAT=-3"  : "Triangular filter, WLINF-WLSUP to WLINF+WLSUP",
                  "ISAT=-2"  : "Flat filter, WLINF-0.5*WLSUP to WLINF+0.5*WLSUP",
                  "ISAT=-1"  : "User defined, read from filter.dat",
                  "ISAT=0"   : "WLINF to WLSUP with filter function = 1 (default)",
                  "ISAT=1"   : "METEO",
                  "ISAT=2"   : "GOES(EAST)",
                  "ISAT=3"   : "GOES(WEST)",
                  "ISAT=4"   : "AVHRR1(NOAA8)",
                  "ISAT=5"   : "AVHRR2(NOAA8)",
                  "ISAT=6"   : "AVHRR1(NOAA9)",
                  "ISAT=7"   : "AVHRR2(NOAA9)",
                  "ISAT=8"   : "AVHRR1(NOAA10)",
                  "ISAT=9"   : "AVHRR2(NOAA10)",
                  "ISAT=10"  : "AVHRR1(NOAA11)",
                  "ISAT=11"  : "AVHRR2(NOAA11)",
                  "ISAT=12"  : "GTR-100 ch1",
                  "ISAT=13"  : "GTR-100 ch2",
                  "ISAT=14"  : "GTR-100 410nm channel",
                  "ISAT=15"  : "GTR-100 936nm channel",
                  "ISAT=16"  : "MFRSR 415nm channel",
                  "ISAT=17"  : "MFRSR 500nm channel",
                  "ISAT=18"  : "MFRSR 610nm channel",
                  "ISAT=19"  : "MFRSR 665nm channel",
                  "ISAT=20"  : "MFRSR 862nm channel",
                  "ISAT=21"  : "MFRSR 940nm channel",
                  "ISAT=22"  : "AVHRR3 (nominal)",
                  "ISAT=23"  : "AVHRR4 (nominal)",
                  "ISAT=24"  : "AVHRR5 (nominal)",
                  "ISAT=25"  : "Biological action spectra for DNA damage by UVB radiation",
                  "ISAT=26"  : "AIRS1 380-460nm",
                  "ISAT=27"  : "AIRS2 520-700nm",
                  "ISAT=28"  : "AIRS3 670-975nm",
                  "ISAT=29"  : "AIRS4 415-1110nm",
                  "NOTHRM=-1": "Thermal emission turned on when lambda > 2um",
                  "NOTHRM=0" : "Thermal emission turned on for all wavelengths",
                  "NOTHRM=1" : "No thermal emission"
                 }
        if key in docstr:
            return "{:10s}# {}".format(key, docstr[key])
        else:
            return key

    def setRtRange(self):
        '''
        rt parameter range dictionary.  key is rt parameter name, value is description and suggested range and skew
        :return:
        '''
        return {"SZA":    "Solar zenith angle (degrees)            $ 0.0:90:1",
                "CSZA":   "Cosine of solar zenith                  $ 0:1:1",
                "WLINF":  "Wavelength lower limit (um)             $ 0.250:100:1",
                "WLSUP":  "Wavelength upper limit (um)             $ 0.250:100:1",
                "WLINC":  "Wavelength/wavenumber increment         $ -.01;0;20",
                "CSZA":   "cosine of solar zenith                  $ 0:1:1",
                "ZPRES":  "Effective surface altitude (km)         $ 0:5:1",
                "PBAR":   "Surface pressure (mb)                   $ 500:2000:1",
                "SCLH2O": "Water vapor scale height (km)           $ 0.5:5:1",
                "UW":     "Integrated water vapor (g/cm2)          $ 0:8:0.3",
                "UO3":    "Integrated ozone amount (atm-cm)        $ 0.050:0.200:0.3",
                "O3TRP":  "Tropospheric ozone amount (atm-cm)      $ 0.00:0.02:0.3",
                "ZTRP":   "Tropospheric altitude (km)              $ 5:12:1",
                "XRSC":   "Rayleigh scattering sensitivity factor  $ 0.5:2:1",
                "XN2":    "N2 volume mixing ratio (ppm)            $ 5e5:1e6:1",
                "XO2":    "O2 volume mixing ratio (ppm)            $ 1e5:4e5:1",
                "XCO2":   "CO2 volume mixing ratio (ppm)           $ 0:800:1",
                "XCH4":   "CH4 volume mixing ratio (ppm)           $ 0:3.48:1",
                "XN2O":   "N2O volume mixing ratio (ppm)           $ 0:0.64:1",
                "XCO":    "CO volume mixing ratio (ppm)            $ 0:30:1",
                "XNO2":   "NO2 volume mixing ratio (ppm)           $ 0:4.6e-5:1",
                "XSO2":   "SO2 volume mixing ratio (ppm)           $ 0:6e-4:1",
                "XNH3":   "NH3 volume mixing ratio (ppm)           $ 0:1.0e-3:1",
                "XNO":    "NO volume mixing ratio (ppm)            $ 0:6.0e-4:1",
                "XHNO3":  "HNO3 volume mixing ratio (ppm)          $ 0:1e-6:1",
                "XO4":    "O4 density sensitivity factor           $ 0.5:2.0:1",
                "ALBCON": "Surface albedo                          $ 0:1:1",
                "ISALB":  "Surface albedo model                    $ -1;0;1;2;3;4;5;6;7;8;9;10",
                "ISAT":   "Filter function                         $ -4;-3;-2;-1;0;1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29",
                "ZCLOUD": "Cloud height (km,km,km,km,km)           $ 0:10:1",
                "TCLOUD": "Cloud optical depth (,,,,,)             $ 0:100:0.3",
                "LWP":    "Liquid water path (g/m2)                $ 0:1000:0.3",
                "NRE":    "Cloud drop radius (um)                  $ 2:128:0",
                "RHCLD":  "Cloud relative humidity                 $ 0.0:1.0:1",
                "JAER":   "Stratospheric Aerosol type              $ 0;1;2;3;4",
                "ZAER":   "SA layer altitudes (km,km,km,km,km)     $ 0:10:1",
                "TAERST": "SA Optical depths (,,,,,)               $ 0:10:0.3",
                "IAER":   "Bounday Layer aerosol type              $ 0;1;2;3;4",
                "RHAER":  "Relative humidity for BL aerosols       $ 0:1:1",
                "TBAER":  "Optical depth of BL aerosols            $ 0:10:0.3",
                "BTEMP":  "Surface termperature (K)                $ 180:330:1",
                "NOTHRM": "Thermal calculation selector            $ -1;0;1",
                "IOUT":   "Output format selector                  $ 1;2;10;11;20;21",
                "ZOUT":   "Output altitude (km,km)                 $ 0.0:10.0:1"}

