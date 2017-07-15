import numpy as np
def NiceStep(self, top):
    v=np.log10(abs(top))
    mantissa = int(v)
    fraction = v-mantissa
    step = 1
    if fraction < 0.30103:
        step = 0.1
    elif fraction < 0.69897:
        step = 0.2
    else:
        step = 0.5
    return step*10**mantissa

NiceStep(0,10)


