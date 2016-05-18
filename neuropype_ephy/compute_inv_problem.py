# 2015.10.09 12:02:45 EDT
# Embedded file name: /home/karim/Documents/pasca/packages/neuropype_ephy/neuropype_ephy/compute_inv_problem.py
"""
Created on Thu Oct  8 17:53:07 2015

@author: pasca
"""


# compute noise covariance data from a continuous segment of raw data.
# Employ empty room data (collected without the subject) to calculate
# the full noise covariance matrix.
# This is recommended for analyzing ongoing spontaneous activity.
def compute_noise_cov(cov_fname, raw):
    import os.path as op

    from mne import compute_raw_covariance, pick_types, write_cov
    from nipype.utils.filemanip import split_filename as split_f

    print '***** COMPUTE RAW COV *****' + cov_fname

    if not op.isfile(cov_fname):

        data_path, basename, ext = split_f(raw.info['filename'])
        fname = op.join(data_path, '%s-cov.fif' % basename)

        reject = dict(mag=4e-12, grad=4000e-13, eog=250e-6)

        picks = pick_types(raw.info, meg=True, ref_meg=False, exclude='bads')

        noise_cov = compute_raw_covariance(raw, picks=picks, reject=reject)

        write_cov(fname, noise_cov)

    else:
        print '*** NOISE cov file %s exists!!!' % cov_fname

    return cov_fname


def read_noise_cov(cov_fname, raw_info):
    import os.path as op
    import numpy as np
    import mne

    print '***** READ RAW COV *****' + cov_fname

    if not op.isfile(cov_fname):
        # create an Identity matrix
        picks = mne.pick_types(raw_info, meg=True, ref_meg=False,
                               exclude='bads')
        ch_names = [raw_info['ch_names'][i] for i in picks]

        C = mne.Covariance(data=np.identity(len(picks)), names=ch_names,
                           bads=[], projs=[], nfree=0)
        mne.write_cov(cov_fname, C)
    else:
        print '*** noise covariance file %s exists!!!' % cov_fname
        noise_cov = mne.read_cov(cov_fname)

    return noise_cov


# compute inverse solution on raw data
def compute_ts_inv_sol(raw, fwd_filename, cov_fname, snr, inv_method, aseg):
    import os.path as op
    import numpy as np
    import mne
    from mne.minimum_norm import make_inverse_operator, apply_inverse_raw
    from nipype.utils.filemanip import split_filename as split_f

    print '***** READ FWD SOL %s *****' % fwd_filename
    forward = mne.read_forward_solution(fwd_filename)

    # Convert to surface orientation for cortically constrained
    # inverse modeling
    if not aseg:
        forward = mne.convert_forward_solution(forward, surf_ori=True)

    lambda2 = 1.0 / snr ** 2

    # compute inverse operator
    print '***** COMPUTE INV OP *****'
    inverse_operator = make_inverse_operator(raw.info, forward, cov_fname,
                                             loose=0.2, depth=0.8)

    # apply inverse operator to the time windows [t_start, t_stop]s
    # TEST
    t_start = 0  # sec
    t_stop = 3  # sec
    start, stop = raw.time_as_index([t_start, t_stop])
    print '***** APPLY INV OP ***** [%d %d]sec' % (t_start, t_stop)
    stc = apply_inverse_raw(raw, inverse_operator, lambda2, inv_method,
                            label=None,
                            start=start, stop=stop, pick_ori=None)

    print '***'
    print 'stc dim ' + str(stc.shape)
    print '***'

    subj_path, basename, ext = split_f(raw.info['filename'])
    data = stc.data

    print 'data dim ' + str(data.shape)

    # save results in .npy file that will be the input for spectral node
    print '***** SAVE SOL *****'
    ts_file = op.abspath(basename + '.npy')
    np.save(ts_file, data)

    return ts_file

'''
+---------------------+-----------+-----------+-----------+-----------------+--------------+
| Inverse desired                             | Forward parameters allowed                 |
+=====================+===========+===========+===========+=================+==============+
|                     | **loose** | **depth** | **fixed** | **force_fixed** | **surf_ori** |
+---------------------+-----------+-----------+-----------+-----------------+--------------+
| | Loose constraint, | 0.2       | 0.8       | False     | False           | True         |
| | Depth weighted    |           |           |           |                 |              |
+---------------------+-----------+-----------+-----------+-----------------+--------------+
| | Loose constraint  | 0.2       | None      | False     | False           | True         |
+---------------------+-----------+-----------+-----------+-----------------+--------------+
| | Free orientation, | None      | 0.8       | False     | False           | True         |
| | Depth weighted    |           |           |           |                 |              |
+---------------------+-----------+-----------+-----------+-----------------+--------------+
| | Free orientation  | None      | None      | False     | False           | True | False |
+---------------------+-----------+-----------+-----------+-----------------+--------------+
| | Fixed constraint, | None      | 0.8       | True      | False           | True         |
| | Depth weighted    |           |           |           |                 |              |
+---------------------+-----------+-----------+-----------+-----------------+--------------+
| | Fixed constraint  | None      | None      | True      | True            | True         |
+---------------------+-----------+-----------+-----------+-----------------+--------------+       
'''


# compute the inverse solution on raw data considering N_r regions in source
# space  based on a FreeSurfer cortical parcellation
def compute_ROIs_inv_sol(raw_filename, sbj_id, sbj_dir, fwd_filename, cov_fname, snr,
                         inv_method, parc, aseg, aseg_labels):
    import os.path as op
    import numpy as np
    import mne
    import pickle

    from mne.io import Raw
    from mne.minimum_norm import make_inverse_operator, apply_inverse_raw
    from nipype.utils.filemanip import split_filename as split_f

    from neuropype_ephy.compute_inv_problem import get_aseg_labels

    print '***** READ raw filename %s *****' % raw_filename
    raw = Raw(raw_filename)

    print '***** READ noise covariance %s *****' % cov_fname
    noise_cov = mne.read_cov(cov_fname)

    print '***** READ FWD SOL %s *****' % fwd_filename
    forward = mne.read_forward_solution(fwd_filename)

    print '***** SNR %s *****' % snr
    
    if not aseg:
        forward = mne.convert_forward_solution(forward, surf_ori=True,
                                               force_fixed=False)

    lambda2 = 1.0 / snr ** 2

    # compute inverse operator
    print '***** COMPUTE INV OP *****'
    if not aseg:
        loose = 0.2
        depth = 0.8
    else:
        loose = None
        depth = None

    inverse_operator = make_inverse_operator(raw.info, forward, noise_cov,
                                             loose=loose, depth=depth,
                                             fixed=False)

    # apply inverse operator to the time windows [t_start, t_stop]s
    print '***** APPLY INV OP *****'
    stc = apply_inverse_raw(raw, inverse_operator, lambda2, inv_method,
                            label=None,
                            start=None, stop=None,
                            buffer_size=1000,
                            pick_ori=None)  # None 'normal'

    print '***'
    print 'stc dim ' + str(stc.shape)
    print '***'

    labels_cortex = mne.read_labels_from_annot(sbj_id, parc=parc,
                                               subjects_dir=sbj_dir)

    src = inverse_operator['src']

    # allow_empty : bool -> Instead of emitting an error, return all-zero time
    # courses for labels that do not have any vertices in the source estimate
    # TODO cosa accade se la uso con solo la cortex? -> OK!!!
    label_ts = mne.extract_label_time_course_AP(stc, labels_cortex, src,
                                                mode='mean_flip',
                                                allow_empty=True,
                                                return_generator=False)

    # save results in .npy file that will be the input for spectral node
    print '***** SAVE SOL *****'
    subj_path, basename, ext = split_f(raw.info['filename'])
    ts_file = op.abspath(basename + '_ROI_ts.npy')
    np.save(ts_file, label_ts)

    if aseg:
        labels_aseg = get_aseg_labels(src, sbj_dir, sbj_id)
        labels = labels_cortex + labels_aseg
    else:
        labels = labels_cortex

    print labels[0].pos
    print len(labels)

    labels_file = op.abspath('labels.dat')
    with open(labels_file, "wb") as f:
        pickle.dump(len(labels), f)
        for value in labels:
            pickle.dump(value, f)

    label_names_file = op.abspath('label_names.txt')
    label_coords_file = op.abspath('label_coords.txt')

    label_names = []
    label_coords = []

    for value in labels:
        label_names.append(value.name)
#        label_coords.append(value.pos[0])
        label_coords.append(np.mean(value.pos, axis=0))

    np.savetxt(label_names_file, np.array(label_names, dtype=str),
               fmt="%s")
    np.savetxt(label_coords_file, np.array(label_coords, dtype=float),
               fmt="%f %f %f")

    return ts_file, labels_file, label_names_file, label_coords_file


# return a list of Label objs
def get_aseg_labels(src, sbj_dir, sbj_id):
    import os.path as op
    import numpy as np

    from mne import Label
    from mne import get_volume_labels_from_aseg_AP

    # read the aseg file
    aseg_fname = op.join(sbj_dir, sbj_id, 'mri/aseg.mgz')
    all_labels_aseg = get_volume_labels_from_aseg_AP(aseg_fname)

    # creo una lista di label per aseg
    labels_aseg = list()
    for nr in range(2, len(src)):
        vertices = src[nr]['vertno']

        pos = src[nr]['rr'][src[nr]['vertno'], :]
        roi_str = src[nr]['seg_name']
        try:
            ind = all_labels_aseg[0].index(roi_str)
            color = np.array(all_labels_aseg[1][ind])/255
        except ValueError:
            pass

        if 'left' in roi_str.lower():
            hemi = 'lh'
            roi_str = roi_str.replace('Left-', '') + '-lh'
        elif 'right' in roi_str.lower():
            hemi = 'rh'
            roi_str = roi_str.replace('Right-', '') + '-rh'
        else:
            hemi = 'both'

        label = Label(vertices=vertices, pos=pos, hemi=hemi,
                      name=roi_str, color=color,
                      subject=sbj_id)
        labels_aseg.append(label)

    return labels_aseg
