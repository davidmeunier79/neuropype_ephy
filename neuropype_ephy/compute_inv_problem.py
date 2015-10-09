# 2015.10.09 12:02:45 EDT
# Embedded file name: /home/karim/Documents/pasca/packages/neuropype_ephy/neuropype_ephy/compute_inv_problem.py
"""
Created on Thu Oct  8 17:53:07 2015

@author: pasca
"""
# compute inverse solution on raw data
def compute_inv_sol(raw, forward, snr, method):
    import os.path as op
    import numpy as np
    import mne
    from mne.minimum_norm import make_inverse_operator, apply_inverse_raw
    from nipype.utils.filemanip import split_filename as split_f
    
    lambda2 = 1.0 / snr ** 2
    reject = dict(mag=4e-12, grad=4e-10, eog=0.00025)
    
    picks = mne.pick_types(raw.info, meg=True, ref_meg=False, exclude='bads')

    # compute noise covariance data
    noise_cov = mne.compute_raw_data_covariance(raw, picks=picks, reject=reject)
    
    # compute inverse operator    
    inverse_operator = make_inverse_operator(raw.info, forward, noise_cov, loose=0.2, depth=0.8)
    
    # apply inverse operator to the time windows [t_start, t_stop]s
    t_start = 0 # sec
    t_stop  = 3 # sec
    start, stop = raw.time_as_index([t_start, t_stop])
    stc = apply_inverse_raw(raw, inverse_operator, lambda2, method, label=None, start=start, stop=stop, pick_ori=None)
    
    print '***'
    print stc.shape
    print '***'
    
    subj_path, basename, ext = split_f(raw.info['filename'])
    data = stc.data
    
    print data.shape

    # save results in .npy file that will be the input for spectral node
    ts_file = op.abspath(basename + '.npy')
    np.save(ts_file, data)

    return ts_file