import numpy as np

class RtReader:
    """class to parse sbdart output formats"""
    def __init__(self, cmds, boolvar):
        """
        initialize sbdart output parser appropriate for specified value of IOUT (found in cmd)
        :param cmds: a string variable containing sbdart run commands
        :param boolvar: a boolean object used in Tkinter widgets
        """
        self.rtkeys=[]
        self.rtunits={}     # dictionary of rt quantities and there units
        self.rtdefault=""
        self.intensity=False
        self.plottype=""
        self.phi=[]
        self.zen=[]
        self.IOUT=10

        lines = cmds.split('\n')
        for line in lines:
            line = line.split('#')[0].strip()
            if len(line) == 0 or not '=' in line:
                continue
            p,v=line.split('=')
            if p.strip().upper() == 'IOUT':
                self.IOUT=int(v)
        if self.IOUT == 1:
            self.RT_01_setup()
        elif self.IOUT == 2:
            self.RT_02_setup()
        elif self.IOUT == 11:
            self.RT_11_setup()
        elif self.IOUT in (20,21):
            self.RT_20_setup()
        else:
            self.RT_10_setup()
        self.menucheck={}

        for k in self.rtunits.keys():
            self.menucheck[k]=boolvar()
            if self.rtdefault == k:
                self.menucheck[k].set(True)

        #if self.intensity:
        #    self.menucheck[self.intensitykey]=boolvar()

        self.radiovalue=""

    def Plotables(self):
        """
        :return: list of plotable rt quantities
        """
        fkeys = []
        if self.radiovalue:
            fkeys = [self.radiovalue]
        else:
            fkeys = [k for k in self.menucheck.keys() if self.menucheck[k].get() and not k.startswith(" ")]
        return fkeys


    def RT(self, out, key, ydict):

        if self.IOUT == 1:
            self.RT_01(out, key, ydict)
        elif self.IOUT == 2:
            self.RT_02(out, key, ydict)
        elif self.IOUT == 11:
            self.RT_11(out, key, ydict)
        elif self.IOUT in (20,21):
            self.RT_20(out, key, ydict)
        else:
            self.RT_10(out, key, ydict)

    # def IntensityParm(self, parm):
    #     intensify=self.menucheck.has_key(self.intensitykey) and self.menucheck[self.intensitykey].get()
    #     has_m2="m^2" in self.rtunits[parm]
    #     return intensify and has_m2

    def IntensityLabel(self, label):
            return label.replace("m^2", "m^2/\mu m")

    def keymaker(self, fkey, pkey):
        '''
        generate full dictionary key for yvariable
        :param fkey:
        :param pkey:
        :return: fullkey
        '''
        if len(pkey) > 0:
            fullkey = "{} {}".format(fkey, pkey)
        else:
            fullkey = fkey
        return fullkey


    def RT_01_setup(self):
        self.rtkeys = ["WL", "FFEW", "TOPDN", "TOPUP", "TOPDIR", "BOTDN", "BOTUP", "BOTDIR", "TOPFLUX", "BOTFLUX", "ABSORPTION"]
        self.rtunits = {" Wavenumber":"", " EffectiveTemp":"", " Transmitted":"",
                        "TOPDN":"$W/m^2/\mu m$",   "TOPUP": "$W/m^2/\mu m$",  "TOPDIR": "$W/m^2/\mu m$",
                        "BOTDN":"$W/m^2/\mu m$",   "BOTUP":"$W/m^2/\mu m$",   "BOTDIR":"$W/m^2/\mu m$",
                        "TOPFLUX":"$W/m^2/\mu m$", "BOTFLUX":"$W/m^2/\mu m$", "ABSORPTION":"$W/m^2/\mu m$"}
        self.intensity=False
        self.rtdefault="TOPUP"

    def RT_01(self, out, label, ydict):
        """
        spectral output
        set up keys for IOUT=10 output parameters
        :param fout:
        :param label:
        :return:
        """
        fout=[]
        nwl = int(out.split()[1])
        vout = [float(v) for v in out.split()[2:]]
        for k in range(0,8):
            fout.append(vout[k::8])
        fout.append([fout[5][k]-fout[6][k] for k in range(0,nwl)]) # botflx 8
        fout.append([fout[2][k]-fout[3][k] for k in range(0,nwl)]) # topflx 9
        fout.append([fout[9][k]-fout[8][k] for k in range(0,nwl)]) # abs 10

        j = -1
        for fkey in self.rtkeys:
            j += 1
            fullkey = self.keymaker(fkey, label)
            if not fullkey in ydict.keys():
                ydict[fullkey] = fout[j]

    def RT_02_setup(self):
        self.rtkeys = ["WL", "H2O", "CO2", "O3", "N2O", "CO", "CH4", "O2+N2", "TRACE", "TOTAL"]
        self.rtunits = {" Wavenumber":"", "H2O":"", "CO2":"", "O3":"", "N2O":"", "CO":"",
                        "CH4":"", "O2+N2":"", "TRACE":"", "TOTAL":""}
        self.intensity=False
        self.rtdefault="TOTAL"

    def RT_02(self, out, label, ydict):
        """
        spectral output
        set up keys for IOUT=10 output parameters
        :param fout:
        :param label:
        :return:
        """
        fout=[]
        nwl = int(out.split()[1])
        vout = [float(v) for v in out.split()[2:]]
        for k in range(0,10):
            fout.append(vout[k::10])

        j = -1
        for fkey in self.rtkeys:
            j += 1
            fullkey = self.keymaker(fkey, label)
            if not fullkey in ydict.keys():
                ydict[fullkey] = fout[j]


    def RT_10_setup(self):
        self.rtkeys = ["WLINF","WLSUP","FFEW","TOPDN", "TOPUP", "TOPDIR", "BOTDN", "BOTUP", "BOTDIR", "TOPFLUX", "BOTFLUX", "ABSORPTION"]
        self.rtunits = {" Intensity":True, " Transmitted":False,
                        "TOPDN": "$W/m^2$", "TOPUP": "$W/m^2$", "TOPDIR": "$W/m^2$",
                        "BOTDN":"$W/m^2$",   "BOTUP":"$W/m^2$",   "BOTDIR":"$W/m^2$",
                        "TOPFLUX":"$W/m^2$", "BOTFLUX":"$W/m^2$", "ABSORPTION":"$W/m^2$"}
        self.intensity=True
        self.rtdefault="BOTDN"

    def RT_10(self, out, label, ydict):
        """
        spectrally integrated flux
        set up keys for IOUT=10 output parameters
        :param fout:
        :param label:
        :return:
        """
        vout = [float(v) for v in out.split()]
        fout = vout[0:]
        fout.append(vout[3] - vout[4])
        fout.append(vout[6] - vout[7])
        fout.append(vout[3] - vout[4] - vout[6] + vout[7])
        j = -1
        for fkey in self.rtkeys:
            fullkey = self.keymaker(fkey, label)
            if not fullkey in ydict.keys():
                ydict[fullkey] = []
            j += 1
            ydict[fullkey].append(fout[j])


    def RT_11_setup(self):
        self.rtkeys = ["FFEW","ZZ","PP","FXDN","FXUP", "FXDIR", "DFDZ", "HEAT"]
        self.rtunits = {" Intensity":True, "FXDN": "$W/m^2$", "FXUP": "$W/m^2$", "FXDIR": "$W/m^2$","PP":"mb", "DFDZ":"$mW/m^3$", "HEAT":"K/day"}
        self.intensity=False
        self.rtdefault="FXDN"

    def RT_11(self, out, label, ydict):
        """
        vertical profiles
        set up keys for IOUT=10 output parameters
        :param fout:
        :param label:
        :return:
        """
        vout = [float(v) for v in out.split()]
        nz = int(vout[0])
        phidw = float(vout[1])
        fout=[]
        fout.append([phidw]*nz)
        for k in range(2,9):
            fout.append(vout[k::7])

        j = -1
        for fkey in self.rtkeys:
            fullkey = self.keymaker(fkey, label)
            if not fullkey in ydict.keys():
                ydict[fullkey] = []
            j += 1
            ydict[fullkey]=fout[j]

    def RT_20_setup(self):
        self.rtkeys = ["WLINF","WLSUP","FFEW","TOPDN", "TOPUP", "TOPDIR", "BOTDN", "BOTUP", "BOTDIR", "TOPFLUX", "BOTFLUX", "ABSORPTION", "RADIANCE"]
        self.rtunits = {" Intensity":True, "TOPDN": "$str^{-1}$", "TOPUP": "$str^{-1}$", "TOPDIR": "$str^{-1}$", "BOTDN": "$str^{-1}$", "BOTUP": "$str^{-1}$",
                         "BOTDIR":"$str^{-1}$", "TOPFLUX":"$str^{-1}$", "BOTFLUX":"$str^{-1}$", "ABSORPTION":"$str^{-1}$","RADIANCE":"$W/m^2/str$"}
        self.intensity=True
        self.rtdefault="RADIANCE"

    def RT_20(self, out, label, ydict):
        """
        radiance output
        set up keys for IOUT=20 output parameters
        :param fout:
        :param label:
        :return:
        """
        vout = [float(v) for v in out.split()]
        fout = vout[0:9]
        fout.append(vout[3] - vout[4])
        fout.append(vout[6] - vout[7])
        fout.append(vout[3] - vout[4] - vout[6] + vout[7])

        nphi=int(vout[9])
        nzen=int(vout[10])
        azm1=[float(x) for x in vout[11:11+nphi]]
        azm2=[360-x for x in azm1[-2::-1]]
        self.phi = azm1 + azm2
        self.zen=[float(x) for x in vout[11+nphi:11+nphi+nzen]]
        #print "phi: "," ".join([str(x) for x in self.phi])
        #print "zen: "," ".join([str(x) for x in self.zen])
        nr=11+nphi+nzen
        # mirror image to fill in 180 to 360
        radv=[]

        for k in range(0,nzen):
            i=nr+k*nphi
            r=vout[i:i+nphi]
            radv+=r+r[-2::-1]
        rad=np.array(radv)
        rad.shape = (nzen,2*nphi-1)   # rad[i,j] => ith phi, jth zenith
        rad = np.transpose(rad)
        #print "rad:", type(rad)
        #self.PrintArray(rad, width=8)
        fout.append(rad)
        j = -1
        for fkey in self.rtkeys:
            fullkey = self.keymaker(fkey, label)
            if not fullkey in ydict.keys():
                ydict[fullkey] = []
            j += 1
            ydict[fullkey]=fout[j]

    def PrintArray(self,array,**kwargs):
        width=10
        if kwargs.has_key('width'):
            width=kwargs['width']
        ifmt = "{:"+str(width)+"}"
        gfmt = "{:"+str(width)+".3g}"
        nx,ny=array.shape
        line="".join(ifmt.format(x) for x in range(1,nx+1))
        print line
        for iy in range(0,ny):
            iys="{:3d}".format(iy+1)
            line=iys+"".join([gfmt.format(x) for x in array[:,iy]])
            print line
