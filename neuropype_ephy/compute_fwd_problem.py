# 2015.10.09 12:03:07 EDT
# Embedded file name: /home/karim/Documents/pasca/packages/neuropype_ephy/neuropype_ephy/compute_fwd_problem.py
"""
Created on Mon Oct  5 17:36:56 2015

@author: pasca

Compute leadfield matrix by BEM
"""

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
        fwd_filename = op.join(data_path, '%s-%s-aseg-fwd.fif' % (raw_fname, spacing))  
    else:
        fwd_filename = op.join(data_path, '%s-%s-fwd.fif' % (raw_fname, spacing))
    
    # check if we have just created the fwd matrix    
    if not op.isfile(fwd_filename):            
        # check if bem-sol was created, if not it creates the bem sol using C MNE
        bem_fname   = op.join(bem_dir, '%s-5120-bem-sol.fif' % sbj_id)
        model_fname = op.join(bem_dir, '%s-5120-bem.fif' % sbj_id)   
        if not op.isfile(bem_fname):
            # chek if inner_skull surf exists, if not BEM computation is 
            # performed by MNE python functions mne.bem.make_watershed_bem  
            if not (op.isfile(sbj_inner_skull_fname) or op.isfile(inner_skull_fname)):
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
                                         spacing=spacing.replace('-',''),
                                         add_dist=False, overwrite=True, n_jobs=2)
            print '*** source space file %s written ***' % src_fname                              
        else:
            print '*** source space file %s exists!!!' % src_fname
            src = mne.read_source_spaces(src_fname)
        
        if aseg:
            src = create_mixed_source_space(sbj_dir, sbj_id, spacing, labels, src)
#            src_aseg_fname = op.join(bem_dir, '%s-%s-aseg-src.fif' % (sbj_id, spacing))            
#            aseg_fname = op.join(sbj_dir,sbj_id, 'mri/aseg.mgz')
#            
#            if spacing == 'oct-6':
#                pos = 5.0
#            elif spacing == 'ico-5':
#                pos = 3.0  
#                
#            for l in labels:
#                print l
#                vol_label = mne.setup_volume_source_space(sbj_id, mri=aseg_fname, 
#                                          pos=pos,
#                                          bem=model_fname,
#                                          volume_label=l,
#                                          subjects_dir=sbj_dir)
#                src += vol_label
#
#
#            mne.write_source_spaces(src_aseg_fname, src)
#
#            # Export source positions to nift file
#            nii_fname = op.join(bem_dir, '%s-%s-aseg-src.nii' % (sbj_id, spacing)) 
#        
#            # Combine the source spaces
#            src.export_volume(nii_fname, mri_resolution=True)
        
        n = sum(src[i]['nuse'] for i in range(len(src)))
        print('il src space contiene %d spaces e %d vertici' % (len(src),n))  
       
        # check if the co-registration file was created, if not raise an runtime error    
        trans_fname = op.join(data_path, '%s-trans.fif' % raw_fname)
        if not op.isfile(trans_fname):
            raise RuntimeError('coregistration file %s NOT found!!!' % trans_fname)
            
       
        # if all is ok creates the fwd matrix
        mne.make_forward_solution(raw_info, trans_fname, src, bem, fwd_filename,
                                  mindist=5.0, # ignore sources <= 0mm from inner skull
                                  meg=True, eeg=False,
                                  n_jobs=2,
                                  overwrite=True)

    else:
        print '*** FWD file %s exists!!!' % fwd_filename
        
    return fwd_filename

def create_mixed_source_space(sbj_dir, sbj_id, spacing, labels, src):
    import os.path as op
    from mne import setup_volume_source_space,write_source_spaces 
    
    bem_dir = op.join(sbj_dir, sbj_id, 'bem')
    
    src_aseg_fname = op.join(bem_dir, '%s-%s-aseg-src.fif' % (sbj_id, spacing))            
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


    write_source_spaces(src_aseg_fname, src)

    # Export source positions to nift file
    nii_fname = op.join(bem_dir, '%s-%s-aseg-src.nii' % (sbj_id, spacing)) 

    # Combine the source spaces
    src.export_volume(nii_fname, mri_resolution=True)
    
    return src


# test function -> TODO eliminare!
def test_compute_LF_matrix():
    import os
    import os.path as op
    import nipype.pipeline.engine as pe
    from nipype.interfaces.mne import WatershedBEM
    import mne
    import mne.io as io
    from mne.minimum_norm import make_inverse_operator, apply_inverse_raw
    from mne.report import Report
    from nipype.utils.filemanip import split_filename as split_f
    main_path = '/home/karim/Documents/pasca/data/resting_state/'
    sbj_id = 'K0002'
    sbj_dir = op.join(main_path, 'FSF')
    bem_dir = op.join(sbj_dir, sbj_id, 'bem')
    surface_dir = op.join(sbj_dir, sbj_id, 'bem/watershed')
    data_dir = op.join(main_path, 'MEG')
    raw_fname = op.join(data_dir, '%s/%s_rest_tsss_mc.fif' % (sbj_id, sbj_id))
    raw = io.Raw(raw_fname, preload=True)
    picks = mne.pick_types(raw.info, meg=True, ref_meg=False, exclude='bads')
    raw.filter(l_freq=0.1, h_freq=300, picks=picks, method='iir', n_jobs=2)
    raw.resample(sfreq=300, npad=0)
    report = Report()
    surfaces = [sbj_id + '_brain_surface',
     sbj_id + '_inner_skull_surface',
     sbj_id + '_outer_skull_surface',
     sbj_id + '_outer_skin_surface']
    new_surfaces = ['brain.surf',
     'inner_skull.surf',
     'outer_skull.surf',
     'outer_skin.surf']
    sbj_inner_skull_fname = op.join(bem_dir, sbj_id + '-' + new_surfaces[1])
    inner_skull_fname = op.join(bem_dir, new_surfaces[1])
    if not (op.isfile(sbj_inner_skull_fname) or op.isfile(inner_skull_fname)):
        bem_IF = WatershedBEM()
        bem_IF.inputs.subject_id = sbj_id
        bem_IF.inputs.subjects_dir = sbj_dir
        bem_IF.inputs.atlas_mode = True
        bem_IF.run()
        for i in range(len(surfaces)):
            os.system('cp %s %s' % (op.join(surface_dir, surfaces[i]), op.join(bem_dir, sbj_id + '-' + new_surfaces[i])))

    else:
        print '*** inner skull surface exists!!!'
    bem = op.join(bem_dir, '%s-5120-bem-sol.fif' % sbj_id)
    if not op.isfile(bem):
        os.system('$MNE_ROOT/bin/mne_setup_forward_model --subject ' + sbj_id + ' --homog --surf --ico 4')
    else:
        print '*** BEM solution file exists!!!'
    src_fname = op.join(bem_dir, '%s-ico-5-src.fif' % sbj_id)
    if not op.isfile(src_fname):
        src = mne.setup_source_space(sbj_id, fname=True, spacing='ico5', subjects_dir=sbj_dir, overwrite=True, n_jobs=2)
    else:
        print '*** source space file exists!!!'
        src = mne.read_source_spaces(src_fname)
    trans_fname = op.join(data_dir, '%s/%s-trans.fif' % (sbj_id, sbj_id))
    data_path, basename, ext = split_f(raw_fname)
    fwd_filename = op.join(data_path, '%s-fwd.fif' % basename)
    forward = mne.make_forward_solution(raw_fname, trans_fname, src, bem, fwd_filename, n_jobs=2, overwrite=True)
    forward = mne.convert_forward_solution(forward, surf_ori=True)
    snr = 1.0
    lambda2 = 1.0 / snr ** 2
    method = 'MNE'
    reject = dict(mag=4e-12, grad=4e-10, eog=0.00025)
    noise_cov = mne.compute_raw_data_covariance(raw, picks=picks, reject=reject)
    inverse_operator = make_inverse_operator(raw.info, forward, noise_cov, loose=0.2, depth=0.8)
    start, stop = raw.time_as_index([0, 30])
    stc = apply_inverse_raw(raw, inverse_operator, lambda2, method, label=None, start=start, stop=stop, pick_ori=None)
    print '***'
    stc.shape
    print '***'
    subj_path, basename, ext = split_f(raw_fname)
    stc_filename = op.join(subj_path, basename)
    stc.save(stc_filename)
    report_filename = op.join(subj_path, basename + '-BEM-report.html')
    print report_filename
    report.save(report_filename, open_browser=False, overwrite=True)
    return

# plot sensitivity map
def plot_forward(fwd, sbj_id, sbj_dir):
    import mne
    import matplotlib.pyplot as plt
    
    leadfield = fwd['sol']['data']
    print 'Leadfield size : %d x %d' % leadfield.shape
    grad_map = mne.sensitivity_map(fwd, ch_type='grad', mode='fixed')
    mag_map  = mne.sensitivity_map(fwd, ch_type='mag', mode='fixed')
    picks_meg = mne.pick_types(fwd['info'], meg=True, eeg=False)
    
    fig, axes = plt.subplots(1, 1, figsize=(10, 8), sharex=True)
    fig.suptitle('Lead field matrix (500 dipoles only)', fontsize=14)
    im = axes.imshow(leadfield[picks_meg, :500], origin='lower', aspect='auto', cmap='RdBu_r')
    axes.set_title('meg'.upper())
    axes.set_xlabel('sources')
    axes.set_ylabel('sensors')
    plt.colorbar(im, ax=axes, cmap='RdBu_r')
    plt.show()
    plt.figure()
    plt.hist([grad_map.data.ravel(), mag_map.data.ravel()], bins=20, label=['Gradiometers', 'Magnetometers'], 
              color=['c', 'b'])
    plt.legend()
    plt.title('Normal orientation sensitivity')
    plt.xlabel('sensitivity')
    plt.ylabel('count')
    plt.show()
    
    grad_map.plot(subject=sbj_id, time_label='Gradiometer sensitivity', subjects_dir=sbj_dir, clim='auto')

# plot bem on MRI
def check_bem(sbj_id, sbj_dir, raw_info, trans_fname, report):
    import os.path as op
    import mne
    from mne.viz import plot_bem, plot_trans
    
    from mayavi import mlab
    import numpy as np
    
    ### plot bem surfs to MRI in the 3 different views
    fig1 = plot_bem(subject=sbj_id, subjects_dir=sbj_dir, orientation='axial', show=False)
    fig2 = plot_bem(subject=sbj_id, subjects_dir=sbj_dir, orientation='sagittal', show=False)
    fig3 = plot_bem(subject=sbj_id, subjects_dir=sbj_dir, orientation='coronal', show=False)
    report.add_figs_to_section([fig1, fig2, fig3], captions=['axial view', 'sagittal view', 'coronal view'], section='BEM surfaces')
    
    ### plot bem surf and source space
    bem_fname = op.join(sbj_dir, sbj_id + '/bem/%s-5120-bem-sol.fif' % sbj_id)
    surf = mne.read_bem_surfaces(bem_fname, patch_stats=True)
    print 'Number of surfaces : %d' % len(surf)
    brain_col = (0.95, 0.83, 0.83)
    cortex_col = (0.91, 0.89, 0.67)
    points = surf[0]['rr']
    faces = surf[0]['tris']
    fig4 = mlab.figure(size=(600, 600), bgcolor=(0, 0, 0))
    mlab.triangular_mesh(points[:, 0], points[:, 1], points[:, 2], faces, color=brain_col, opacity=0.3)
    
    src_fname = op.join(sbj_dir, sbj_id + '/bem/%s-ico-5-src.fif' % sbj_id)
    src = mne.read_source_spaces(src_fname)
    lh_points = src[0]['rr']
    lh_faces = src[0]['tris']
    rh_points = src[1]['rr']
    rh_faces = src[1]['tris']
    mlab.triangular_mesh(lh_points[:, 0], lh_points[:, 1], lh_points[:, 2], lh_faces, color=cortex_col, opacity=0.8)
    mlab.triangular_mesh(rh_points[:, 0], rh_points[:, 1], rh_points[:, 2], rh_faces, color=cortex_col, opacity=0.8)
    picks = mne.pick_types(raw_info, meg=True, ref_meg=False, eeg=False, stim=False, eog=False, exclude='bads')

    ### plot sensors    
    sens_loc = [ raw_info['chs'][picks[i]]['loc'][:3] for i in range(len(picks)) ]
    sens_loc = np.array(sens_loc)
    mlab.points3d(sens_loc[:, 0], sens_loc[:, 1], sens_loc[:, 2], color=(1, 1, 1), opacity=1, scale_factor=0.005)
    
    report.add_figs_to_section(fig4, captions=['source space'], section='BEM cortex sensors')
    fig5 = plot_trans(raw_info, trans_fname, subject=sbj_id, subjects_dir=sbj_dir)
    report.add_figs_to_section(fig5, captions=['MEG <-> MRI coregistration quality'], section='MEG <-> MRI')


if __name__ == '__main__':
    test_compute_LF_matrix()

