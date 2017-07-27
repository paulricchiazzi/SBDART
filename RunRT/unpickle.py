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
    if not file:
        files = os.listdir('RUNS')
        fdict = {}
        for f in files:
            if f.endswith('.pkl'):
                fn = 'RUNS{}{}'.format(os.path.sep, f)
                mtime = os.path.getctime(fn)
                fdict[mtime] = fn
        for i,k in enumerate(sorted(fdict.keys())):
            date = datetime.fromtimestamp(k).strftime('%Y-%m-%d %H:%M:%S')
            print '{:5} {:60} {}'.format(i,fdict[k],date)
        print ''
        answ = raw_input("Enter number or return to cancel: ")

        if answ.isdigit() and int(answ) in range(len(fdict)):
            j = int(answ)
            file = fdict[sorted(fdict.keys())[j]]
        else:
            return [None, None]
    with open(file, 'rb') as fh:
        x,y = pickle.load(fh)
        x=np.array(x).astype(float)
        print y.keys()
        return x,y
