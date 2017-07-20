import numpy as np

class SolarEphemeris:
    '''
    purpose:  compute the solar zenith and azimuth angles and solar flux
              multiplier for a given location, time and day.

    input:
      iday    day of year (used as a fraction of 365.24 days year)

      time    universal time in decimal hours

      alat    geographic latitude of point on earth's surface (degrees)

      alon    geographic longitude of point on earth's surface (degrees)

    output:

      zenith  solar zenith angle (degrees)

      azimuth solar azimuth measured clockwise from due north (degrees)

      solfac  solar flux multiplier. solfac=1./rsun**2
              where rsun is the current earth-sun distance in
              astronomical units.

    example:
        import numpy as np
        import matplotlib.pyplot as plt
        import Zensun
        ephem = Zensun.Zensun(34.0,-134.0)
        h=np.arange(24)
        zen1,phi,solfac=ephem.Zensun(34,h)
        plt.plot(h,zen1)


    '''
    def __init__(self, alat, alon):
        '''
        :param alat:    latitude of observer (deg)
        :param alon:    longitude of observer (deg)
        '''

        self.alat = alat
        self.alon = alon
        self.degpday = 360.0/365.242
        self.eccen = 0.01671
        self.dayph = 2.0              # day of perihelion
        self.dtor = np.pi/180

        self.nday = (
             1,   6,  11,  16,  21,  26,  31,  36,  41,
             46,  51,  56,  61,  66,  71,  76,  81,  86,
             91,  96, 101, 106, 111, 116, 121, 126, 131,
             136, 141, 146, 151, 156, 161, 166, 171, 176,
             181, 186, 191, 196, 201, 206, 211, 216, 221,
             226, 231, 236, 241, 246, 251, 256, 261, 266,
             271, 276, 281, 286, 291, 296, 301, 306, 311,
             316, 321, 326, 331, 336, 341, 346, 351, 356,
             361, 366 )

        self.eqt = (
            -3.23, -5.49, -7.60, -9.48,-11.09,-12.39,-13.34,-13.95,
            -14.23,-14.19,-13.85,-13.22,-12.35,-11.26,-10.01, -8.64,
            -7.18, -5.67, -4.16, -2.69, -1.29, -0.02,  1.10,  2.05,
            2.80,  3.33,  3.63,  3.68,  3.49,  3.09,  2.48,  1.71,
            0.79, -0.24, -1.33, -2.41, -3.45, -4.39, -5.20, -5.84,
            -6.28, -6.49, -6.44, -6.15, -5.60, -4.82, -3.81, -2.60,
            -1.19,  0.36,  2.03,  3.76,  5.54,  7.31,  9.04, 10.69,
            12.20, 13.53, 14.65, 15.52, 16.12, 16.41, 16.36, 15.95,
            15.19, 14.09, 12.67, 10.93,  8.93,  6.70,  4.32,  1.86,
            -0.62, -3.23)

        self.dec = (
             -23.06,-22.57,-21.91,-21.06,-20.05,-18.88,-17.57,-16.13,
             -14.57,-12.91,-11.16, -9.34, -7.46, -5.54, -3.59, -1.62,
             0.36,  2.33,  4.28,  6.19,  8.06,  9.88, 11.62, 13.29,
             14.87, 16.34, 17.70, 18.94, 20.04, 21.00, 21.81, 22.47,
             22.95, 23.28, 23.43, 23.40, 23.21, 22.85, 22.32, 21.63,
             20.79, 19.80, 18.67, 17.42, 16.05, 14.57, 13.00, 11.33,
             9.60,  7.80,  5.95,  4.06,  2.13,  0.19, -1.75, -3.69,
             -5.62, -7.51, -9.36,-11.16,-12.88,-14.53,-16.07,-17.50,
             -18.81,-19.98,-20.99,-21.85,-22.52,-23.02,-23.33,-23.44,
             -23.35,-23.06)

    def suntimes(self, iday):
        '''
        get times of local noon and day duration
        :param iday:
        :return:
         noon       universal time at which sun is highest in the sky
         daylight   total number of daylight hours
        '''
        day, sunlat, sunlon = self.solarcoorinates(iday, 0)

        t0=(90.-self.alat)*self.dtor
        t1=(90.-sunlat)*self.dtor
        p0=self.alon*self.dtor

        #
        # noon occurs when sin(t0)*sin(t1)*cos(p1-p0) is maximized
        #
        # noon is when p1=p0       for sin(t0)*sin(t1)>0
        #           or p1=p0+pi    for sin(t0)*sin(t1)<0


        # sunset occurs when 0=cos(t0)*cos(t1)+sin(t0)*sin(t1)*cos(p1-p0)
        # cos(p1-p0)=-cos(t0)*cos(t1)/(sin(t0)*sin(t1)
        # p1-p0 = acos(1/(tan(t0)*tan(t1))
        # duration =
        # p1=(-15. * (time - 12. + eqtime / 60.))*pi/180
        # time =12-p1*180/(15*pi)-eqtime/60

        sinfac = np.sin(t0)*np.sin(t1)
        cosfac = np.cos(t0)*np.cos(t1)  # cosfac>0  => summer in observer hemisphere
        if sinfac > 0:
            p1=p0
        else:
            p1=p0+np.pi if p0 < np.pi else p0-np.pi
        noon = (sunlon-p1/self.dtor)/15

        tanfac = 1.0/(np.tan(t0)*np.tan(t1))
        if tanfac > 1:
            daylight = 24
        elif tanfac < -1:
            daylight = 0
        else:
            pdiff = 2*abs(np.arccos(tanfac))
            if (pdiff < np.pi) == (cosfac > 0):   # daylight hours greater than 12 in summer hemisphere
                pdiff = 2*np.pi - pdiff
            daylight = (pdiff*180/np.pi)*24/360
        return noon, daylight

    def sunpos(self, iday, time):
        '''
        compute solar ephemeris
        :param iday:    day of year
        :param time:    time of day (could be numpy array)
        :return:
        zenith      solar zenith (radians)
        azimuth     solar azimuth (radians)
        solfac      solar irradiance correction factor for distance
        '''
        day, sunlat, sunlon = self.solarcoorinates(iday, time)

        t0=(90.-self.alat)*self.dtor
        t1=(90.-sunlat)*self.dtor
        p0=self.alon*self.dtor
        p1=sunlon*self.dtor
        zz=np.cos(t0)*np.cos(t1)+np.sin(t0)*np.sin(t1)*np.cos(p1-p0)
        xx=np.sin(t1)*np.sin(p1-p0)
        yy=np.sin(t0)*np.cos(t1)-np.cos(t0)*np.sin(t1)*np.cos(p1-p0)
        azimuth=np.arctan2(xx,yy)
        zenith=np.arccos(zz)
        rsun=1.-self.eccen*np.cos(self.degpday*(day-self.dayph)*self.dtor)
        solfac=1/rsun**2
        return zenith, azimuth, solfac

    def solarcoorinates(self, iday, time):
        day = ((iday - 1) % 365) + 1
        for i, d in enumerate(self.nday):
            if d > iday: break
        frac = (day - self.nday[i-1]) / (self.nday[i] - self.nday[i-1])
        eqtime = self.eqt[i-1] * (1. - frac) + frac * self.eqt[i]
        decang = self.dec[i-1] * (1. - frac) + frac * self.dec[i]
        sunlat = decang
        sunlon = -15. * (time - 12. + eqtime / 60.)
        return day, sunlat, sunlon

if __name__ == "__main__":
    import numpy as np
    import datetime
    import matplotlib.pyplot as plt
    ephem = SolarEphemeris(55.0, 0)
    iday = 337
    noon, daylight = ephem.suntimes(iday)
    h=np.linspace(noon-0.5*daylight, noon+0.5*daylight,30)
    zen,phi,solfac=ephem.sunpos(iday, h)

    # for d in range(1,365,7):
    #     date = datetime.date.fromordinal(d)
    #     noon, daylight = ephem.suntimes(d)
    #     print 'day={:<10}   date={}   noon={:.2f}   daylight={:.2f}'.format(d,date, noon,daylight)
    #
    date = datetime.date.fromordinal(iday)
    date = str(date)[5:]
    plt.plot(h,solfac)
    plt.xlabel('UTC hours')
    plt.ylabel('Irradiance Factor')
    plt.title('{:.2f} daylight hours for day {} ({})'.format(daylight, iday, date))
    plt.show()
