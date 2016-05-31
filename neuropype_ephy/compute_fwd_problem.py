# 2015.10.09 12:03:07 EDT
# Embedded file name: /home/karim/Documents/pasca/packages/neuropype_ephy/neuropype_ephy/compute_fwd_problem.py
"""
Created on Mon Oct  5 17:36:56 2015

@author: pasca

Compute leadfield matrix by BEM
"""

# TODO eliminare!!!
# compute LF matrix using mne_setup_forward_model (C MNE)
def compute_LF_matrix(sbj_id, sbj_dir, raw_info, aseg, spacing, labels):
    import os.path as op
    import mne

    from mne.bem import make_watershed_bem
    from mne.report import Report

    from nipype.utils.filemanip import split_filename as split_f

    from neuropype_ephy.compute_fwd_problem import create_mixed_source_space

    report = Report()

    bem_dir = op.join(sbj_dir, sbj_id, 'bem')

    surf_name = 'inner_skull.surf'
    sbj_inner_skull_fname = op.join(bem_dir, sbj_id + '-' + surf_name)
    inner_skull_fname = op.join(bem_dir, surf_name)

    data_path, raw_fname, ext = split_f(raw_info['filename'])

    if aseg:
        fwd_filename = op.join(data_path, '%s-%s-aseg-fwd.fif'
                               % (raw_fname, spacing))
    else:
        fwd_filename = op.join(data_path, '%s-%s-fwd.fif'
                               % (raw_fname, spacing))

    # check if we have just created the fwd matrix
    if not op.isfile(fwd_filename):
        # check if bem-sol was created, if not creates the bem sol using C MNE
        bem_fname = op.join(bem_dir, '%s-5120-bem-sol.fif' % sbj_id)
        model_fname = op.join(bem_dir, '%s-5120-bem.fif' % sbj_id)
        if not op.isfile(bem_fname):
            # chek if inner_skull surf exists, if not BEM computation is
            # performed by MNE python functions mne.bem.make_watershed_bem
            if not (op.isfile(sbj_inner_skull_fname) or
                    op.isfile(inner_skull_fname)):
                print sbj_inner_skull_fname + '---> FILE NOT FOUND!!! ---> BEM is computed'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                make_watershed_bem(sbj_id, sbj_dir, overwrite=True)
            else:
                print '*** inner skull surface exists!!!'

            # Create a BEM model for a subject
            surfaces = mne.make_bem_model(sbj_id, ico=4, conductivity=[0.3],
                                          subjects_dir=sbj_dir)
            # Write BEM surfaces to a fiff file
            mne.write_bem_surfaces(model_fname, surfaces)

            # Create a BEM solution using the linear collocation approach
            bem = mne.make_bem_solution(surfaces)
            mne.write_bem_solution(bem_fname, bem)

            print '*** BEM solution file %s written ***' % bem_fname
            # add BEM figures to a Report
            report.add_bem_to_section(subject=sbj_id, subjects_dir=sbj_dir)
            report_filename = op.join(bem_dir, "BEM_report.html")
            print report_filename
            report.save(report_filename, open_browser=False, overwrite=True)
        else:
            bem = bem_fname
            print '*** BEM solution file %s exists!!!' % bem_fname

        # check if source space exists, if not it creates using mne-python fun
        # we have to create the cortical surface source space even when aseg is
        # True
        src_fname = op.join(bem_dir, '%s-%s-src.fif' % (sbj_id, spacing))
        if not op.isfile(src_fname):
            src = mne.setup_source_space(sbj_id, subjects_dir=sbj_dir,
                                         fname=True,
                                         spacing=spacing.replace('-', ''),
                                         add_dist=False, overwrite=True,
                                         n_jobs=2)
            print '*** source space file %s written ***' % src_fname         
        else:
            print '*** source space file %s exists!!!' % src_fname
            src = mne.read_source_spaces(src_fname)

        if aseg:
            src = create_mixed_source_space(sbj_dir, sbj_id, spacing,
                                            labels, src)

        n = sum(src[i]['nuse'] for i in range(len(src)))
        print('il src space contiene %d spaces e %d vertici' % (len(src), n))

        # check if the co-registration file was created
        # if not raise an runtime error
        trans_fname = op.join(data_path, '%s-trans.fif' % raw_fname)
        if not op.isfile(trans_fname):
            raise RuntimeError('coregistration file %s NOT found!!!'
                               % trans_fname)

        # if all is ok creates the fwd matrix
        mne.make_forward_solution(raw_info, trans_fname, src, bem,
                                  fwd_filename,
                                  mindist=5.0, # ignore sources <= 0mm from inner skull
                                  meg=True, eeg=False,
                                  n_jobs=2,
                                  overwrite=True)

    else:
        print '*** FWD file %s exists!!!' % fwd_filename

    return fwd_filename


def create_bem_sol(sbj_dir, sbj_id):
    import os.path as op
    import mne

    from mne.bem import make_watershed_bem
    from mne.report import Report

    report = Report()

    bem_dir = op.join(sbj_dir, sbj_id, 'bem')

    surf_name = 'inner_skull.surf'
    sbj_inner_skull_fname = op.join(bem_dir, sbj_id + '-' + surf_name)
    inner_skull_fname = op.join(bem_dir, surf_name)

    # check if bem-sol was created, if not creates the bem sol using C MNE
    bem_fname = op.join(bem_dir, '%s-5120-bem-sol.fif' % sbj_id)
    model_fname = op.join(bem_dir, '%s-5120-bem.fif' % sbj_id)

    if not op.isfile(bem_fname):
        # chek if inner_skull surf exists, if not BEM computation is
        # performed by MNE python functions mne.bem.make_watershed_bem
        if not (op.isfile(sbj_inner_skull_fname) or
                op.isfile(inner_skull_fname)):
            print inner_skull_fname + '---> FILE NOT FOUND!!!---> BEM computed'
            make_watershed_bem(sbj_id, sbj_dir, overwrite=True)
        else:
            print '\n*** inner skull %s surface exists!!!\n' % inner_skull_fname

        # Create a BEM model for a subject
        surfaces = mne.make_bem_model(sbj_id, ico=4, conductivity=[0.3],
                                      subjects_dir=sbj_dir)

        # Write BEM surfaces to a fiff file
        mne.write_bem_surfaces(model_fname, surfaces)

        # Create a BEM solution using the linear collocation approach
        bem = mne.make_bem_solution(surfaces)
        mne.write_bem_solution(bem_fname, bem)

        print '\n*** BEM solution file %s written ***\n' % bem_fname

        # add BEM figures to a Report
        report.add_bem_to_section(subject=sbj_id, subjects_dir=sbj_dir)
        report_filename = op.join(bem_dir, "BEM_report.html")
        print '\n*** REPORT file %s written ***\n' % report_filename
        print report_filename
        report.save(report_filename, open_browser=False, overwrite=True)
    else:
        bem = bem_fname
        print '\n*** BEM solution file %s exists!!! ***\n' % bem_fname

    return bem


def create_src_space(sbj_dir, sbj_id, spacing, is_blind):
    import os.path as op
    import mne

    bem_dir = op.join(sbj_dir, sbj_id, 'bem')

    # check if source space exists, if not it creates using mne-python fun
    # we have to create the cortical surface source space even when aseg is
    # True
    if is_blind:
        # if is_blind we have to precomputed the source space sincw we had
        # to remove some labels
        src_fname = op.join(bem_dir, '%s-blind-%s-src.fif' % (sbj_id, spacing))
        if not op.isfile(src_fname):
            raise '\n *** you have to compute the source space blind!!! ***\n'
        else:
            print '\n*** source space file %s exists!!!\n' % src_fname
            src = mne.read_source_spaces(src_fname)
    else:
        src_fname = op.join(bem_dir, '%s-%s-src.fif' % (sbj_id, spacing))
        if not op.isfile(src_fname):
            src = mne.setup_source_space(sbj_id, subjects_dir=sbj_dir,
                                         fname=True,
                                         spacing=spacing.replace('-', ''),
                                         add_dist=False, overwrite=True,
                                         n_jobs=2)
            print '\n*** source space file %s written ***\n' % src_fname
        else:
            print '\n*** source space file %s exists!!!\n' % src_fname
            src = mne.read_source_spaces(src_fname)

    return src


def create_mixed_source_space(sbj_dir, sbj_id, spacing, labels, src):
    import os.path as op
    from mne import setup_volume_source_space

    bem_dir = op.join(sbj_dir, sbj_id, 'bem')

#    src_aseg_fname = op.join(bem_dir, '%s-%s-aseg-src.fif' %(sbj_id, spacing))
    aseg_fname = op.join(sbj_dir, sbj_id, 'mri/aseg.mgz')

    if spacing == 'oct-6':
        pos = 5.0
    elif spacing == 'ico-5':
        pos = 3.0

    model_fname = op.join(bem_dir, '%s-5120-bem.fif' % sbj_id)
    for l in labels:
        print l
        vol_label = setup_volume_source_space(sbj_id, mri=aseg_fname,
                                              pos=pos,
                                              bem=model_fname,
                                              volume_label=l,
                                              subjects_dir=sbj_dir)
        src += vol_label

#    write_source_spaces(src_aseg_fname, src)

    # Export source positions to nift file
    nii_fname = op.join(bem_dir, '%s-%s-aseg-src.nii' % (sbj_id, spacing))

    # Combine the source spaces
    src.export_volume(nii_fname, mri_resolution=True)

    return src


def is_trans(raw_info):
    import os.path as op

    from nipype.utils.filemanip import split_filename as split_f

    data_path, raw_fname, ext = split_f(raw_info['filename'])

    # check if the co-registration file was created
    # if not raise an runtime error
    i_ica = raw_fname.find('-preproc')
    if i_ica != -1:
        raw_fname = raw_fname[:i_ica]

    trans_fname = op.join(data_path, '%s-trans.fif' % raw_fname)
    if not op.isfile(trans_fname):
        raise RuntimeError('*** coregistration file %s NOT found!!!'
                           % trans_fname)

    return trans_fname


def compute_fwd_sol(raw_info, trans_fname, src, bem, fwd_filename):
    import mne

    mne.make_forward_solution(raw_info, trans_fname, src, bem,
                              fwd_filename,
                              mindist=5.0, # ignore sources <= 0mm from inner skull
                              meg=True, eeg=False,
                              n_jobs=2,
                              overwrite=True)

    print '\n*** FWD file %s written!!!\n' % fwd_filename
