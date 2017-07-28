'''

use this script to load a RunRT dataset into an ipython session

For example, start ipython and run these commands,

%pylab
run unpickle
x,y = unpickle()
k=y.keys()[4]
plot(x,y[k])


'''

import pickle
import os
import numpy as np
from datetime import datetime

def unpickle(file=''):
    '''
    load a RunRT pickle file into an ipython session
    :param file   name of pickle file, if not specified a menu of possibilities is offered.
    :return:
        x  -- the x vector
        y  -- a dictionary of y vectors
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

            print '{:5}   {:19}  {:>5}  {}'.format(i, date, size, filename[k])
        print ''
        answ = raw_input("Enter number, <Return>=get most recent, q=abort): ")

        if answ.isdigit() and int(answ) in range(len(filename)):
            j = int(answ)
            file = filename[sorted(filename.keys())[j]]
        elif answ == '':
            j=len(filename)-1
            file = filename[sorted(filename.keys())[j]]
        else:
            return [None, None]


    with open(file, 'rb') as fh:
        return pickle.load(fh)

def getkeys(y, pat=''):
    kys = []
    for k in y.keys():
        addit = True
        for p in pat.split('*'):
            if not p in k:
                addit = False
        if addit:
            kys.append(k)
    return kys
