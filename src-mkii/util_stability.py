import numpy as np


def check_unstable(m, cutoff):

    greater = np.abs(m) > cutoff
    if np.any(greater):
        print('MODEL UNSTABLE. EXITING.')
        exit()