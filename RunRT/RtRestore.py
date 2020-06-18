import pickle
import os
import numpy as np
from datetime import datetime
class RtRestore:
    '''
    loads a RunRT dataset into an ipython session

    For example, start ipython and run these commands,

    %pylab
    import RtRestore as u
    p = u.RtRestore()
    ky=p.getkeys()[5]
    p.data[ky]
    '''

    def __init__(self,file=''):
        '''
        load a RunRT pickle file into an ipython session
        :param file   name of pickle file, if not specified a menu of possibilities is offered.

        attributes:
            data  -- a ordered dictionary of y vectors
                  the last key,value pair in the dictionary is the x vector

                  for IOUT = 10, the dependent variable is to the first variant in the RunRT
                            command file. e.g., if TCLOUD is the first variant, its variation
                            is in data['TCLOUD']
                  for IOUT = 1 or 2, the dependent variable is wavelength, given by data['WL']
                  for IOUT = 11 the dependent variable is altitude, given by data['ZZ']
                  for IOUT = 20 or 21 represent output for radiance plots, and do not have a
                            dependent variable forfor line plots
        '''
        if not file:
            files = os.listdir('RUNS')
            filename = {}
            filesize = {}
            for f in files:
                if f.endswith('.pkl'):
                    fn = 'RUNS{}{}'.format(os.path.sep, f)
                    mtime = os.path.getctime(fn)
                    filename[mtime] = fn
                    filesize[mtime] = os.path.getsize(fn)
            for i,k in enumerate(sorted(filename.keys())):
                date = datetime.fromtimestamp(k).strftime('%Y-%m-%d %H:%M:%S')
                size = filesize[k]
                if size < 1000000:
                    size = '{}KB'.format(size/1000)
                elif size < 1000000000:
                    size = '{}MB'.format(size/1000000)
                else:
                    size = '{}GB'.format(size/1000000000)

                print('{:5}   {:19}  {:>5}  {}'.format(i, date, size, filename[k]))
            print('')
            answ = raw_input("Enter number, <Return>=get most recent, q=abort): ")

            if answ.isdigit() and int(answ) in range(len(filename)):
                j = int(answ)
                file = filename[sorted(filename.keys())[j]]
            elif answ == '':
                j=len(filename)-1
                file = filename[sorted(filename.keys())[j]]
            else:
                self.data = None
        with open(file, 'rb') as fh:
            self.data = pickle.load(fh)

    def getkeys(self, pat=''):
        '''
        filter list data dictionary keys according to pat

        :param
            pat: filter pattern used to match against dictionary keys
                 use an asterisk character '*' as a wildcard.

        :return:
            kys: dictionary keys that match pat
        '''
        if not self.data:
            return None
        kys = []
        for k in self.data.keys():
            addit = True
            for p in pat.split('*'):
                if not p in k:
                    addit = False
            if addit:
                kys.append(k)
        return kys
