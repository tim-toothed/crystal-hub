import numpy as np
import warnings

#@cron
def diagonalisation(matrix, wordy=False):
    matrix = np.round(np.copy(matrix),16)
    w,v = np.linalg.eigh(matrix)
    if round(np.linalg.norm(v[:,0]),8) != 1:
        warnings.warn('Not normalized eigenvectors!')
        print('Performing normalization...\n' if wordy else "", end = "")
        for ixx in range(v.shape[1]):
            v[:,ixx] /= np.linalg.norm(v[:,ixx])
        print('...done\n' if wordy else "", end = "")
    w = np.round(w,16)
    v = np.round(v,16)
    return w,v
