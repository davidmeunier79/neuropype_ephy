# -*- coding: utf-8 -*-

import os
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio

from nipype.interfaces.utility import IdentityInterface

from neuropype_ephy.interfaces.mne.LF_computation import LFComputation
from neuropype_ephy.interfaces.mne.Inverse_solution import NoiseCovariance
from neuropype_ephy.interfaces.mne.Inverse_solution import InverseSolution
from neuropype_ephy.pipelines.preproc_meeg import create_pipeline_preproc_meeg
from neuropype_ephy.pipelines.ts_to_conmat import create_pipeline_time_series_to_spectral_connectivity

from neuropype_graph.pipelines.conmat_to_graph import create_pipeline_conmat_to_graph_density

def create_pipeline_source_reconstruction(main_path, sbj_dir,
                                          pipeline_name='inv_sol_pipeline',
                                          spacing='ico-5',
                                          inv_method='MNE',
                                          parc='aparc',
                                          aseg=False,
                                          aseg_labels=[],
                                          noise_cov_fname=None):

    pipeline = pe.Workflow(name=pipeline_name)
    pipeline.base_dir = main_path

    inputnode = pe.Node(IdentityInterface(fields=['sbj_id', 'raw']),
                        name='inputnode')

    # Lead Field computation Node
    LF_computation = pe.Node(interface=LFComputation(), name='LF_computation')

    LF_computation.inputs.sbj_dir = sbj_dir
    LF_computation.inputs.spacing = spacing
    LF_computation.inputs.aseg = aseg
    if aseg:
        LF_computation.inputs.aseg_labels = aseg_labels

    pipeline.connect(inputnode, 'sbj_id', LF_computation, 'sbj_id')

    pipeline.connect(inputnode, ('raw', get_raw_info),
                     LF_computation, 'raw_info')

    # Noise Covariance Matrix Node
    create_noise_cov = pe.Node(interface=NoiseCovariance(),
                               name="create_noise_cov")

    if noise_cov_fname is not None:
        create_noise_cov.inputs.cov_fname_in = noise_cov_fname

    pipeline.connect(inputnode, 'raw', create_noise_cov, 'raw')

    # Inverse Solution Node
    inv_solution = pe.Node(interface=InverseSolution(), name='inv_solution')

    inv_solution.inputs.sbj_dir = sbj_dir
    inv_solution.inputs.inv_method = inv_method
    inv_solution.inputs.parc = parc
    inv_solution.inputs.aseg = aseg
    if aseg:
        inv_solution.inputs.aseg_labels = aseg_labels

    pipeline.connect(inputnode, 'sbj_id', inv_solution, 'sbj_id')
    pipeline.connect(inputnode, 'raw', inv_solution, 'raw')
    pipeline.connect(LF_computation, 'fwd_filename',
                     inv_solution, 'fwd_filename')
    pipeline.connect(create_noise_cov, 'cov_fname_out',
                     inv_solution, 'cov_filename')

    return pipeline


def get_raw_info(raw):
    return raw.info


def get_freq_band(freq_band_name):

    freq_bands = [[8, 12]]  # [15,29],[60,90]
    freq_band_names = ['alpha']  # "beta",'gamma2'
#    from params import freq_band_names,freq_bands

    if freq_band_name in freq_band_names:
        print freq_band_name
        print freq_band_names.index(freq_band_name)

        return freq_bands[freq_band_names.index(freq_band_name)]
        
if __name__ == '__main__':
    main_path = '/home/karim/Documents/Fanny'
    data_path = main_path
    sbj_dir = os.path.join(main_path, 'FSF')

    subject_ids = ['S01']  # 'S02'
    sessions = ['repos_1']  # 'repos_2'

    noise_cov_fname = os.path.join(main_path, 'Big_Noise-cov.fif')

    mod = True
    aseg = True
    
    if mod:
        radatools_optim = "WS trfr 1"
        
#    aseg_labels = ['Left-Accumbens-area',
#                   'Left-Amygdala',
#                   'Left-Caudate',
#                   'Left-Hippocampus',
#                   'Left-Pallidum',
#                   'Left-Putamen',
#                   'Left-Thalamus-Proper',
#                   'Left-Cerebellum-Cortex',
#                   'Brain-Stem',
#                   'Right-Accumbens-area',
#                   'Right-Amygdala',
#                   'Right-Caudate',
#                   'Right-Hippocampus',
#                   'Right-Pallidum',
#                   'Right-Putamen',
#                   'Right-Thalamus-Proper',
#                   'Right-Cerebellum-Cortex']

    aseg_labels = ['Left-Amygdala',
                   'Left-Hippocampus',
                   'Left-Thalamus-Proper',
                   'Left-Cerebellum-Cortex',
                   'Brain-Stem',
                   'Right-Amygdala',
                   'Right-Hippocampus',
                   'Right-Thalamus-Proper',
                   'Right-Cerebellum-Cortex']

    # create main workflow
    main_workflow = pe.Workflow(name='wf_raw_to_graph_analysis')
    main_workflow.base_dir = main_path

    # info source
    infosource = pe.Node(interface=IdentityInterface(fields=['subject_id',
                                                             'sess_index',
                                                             'freq_band_name']),
                         name="infosource")
    infosource.iterables = [('subject_id', subject_ids),
                            ('sess_index', sessions),
                            ('freq_band_name', ['alpha'])]

    # data source
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id',
                                                             'sess_index'],
                                                   outfields=['raw_file']),
                         name='datasource')
    datasource.inputs.base_directory = data_path
    datasource.inputs.template = '%s/%s%s'
    datasource.inputs.template_args = dict(raw_file=[['subject_id',
                                                     'sess_index', ".ds"]])
#    datasource.inputs.template = 'wf_raw_to_graph_analysis/preproc_meeg/_freq_band_name_alpha_sess_index_%s_subject_id_%s/preproc/*_ica%s' # meditation
#    datasource.inputs.template_args = dict(raw_file=[['sess_index',
#                                                     'subject_id', ".npy"]])
    datasource.inputs.sort_filelist = True

    main_workflow.connect(infosource, 'subject_id', datasource, 'subject_id')
    main_workflow.connect(infosource, 'sess_index', datasource, 'sess_index')

    preproc_workflow = create_pipeline_preproc_meeg(main_path,
                                                    is_sensor_space=False,
                                                    is_set_ICA_components=True,
                                                    n_comp_exclude={'S01': [[0,5,53]]},
                                                    data_type='ds')

    main_workflow.connect(datasource, 'raw_file',
                          preproc_workflow, 'inputnode.raw_file')

    if aseg:
        aseg_labels = aseg_labels
    else:
        aseg_labels = []

    inv_sol_workflow = create_pipeline_source_reconstruction(main_path,
                                                             sbj_dir,
                                                             spacing='ico-5',
                                                             aseg = aseg,
                                                             aseg_labels=aseg_labels,
                                                             noise_cov_fname=noise_cov_fname)

    main_workflow.connect(infosource, 'subject_id',
                          inv_sol_workflow, 'inputnode.sbj_id')

#    main_workflow.connect(datasource, 'raw_file',
#                          inv_sol_workflow, 'inputnode.raw')

    main_workflow.connect(preproc_workflow, 'preproc.out_file',
                          inv_sol_workflow, 'inputnode.raw')

    spectral_workflow = create_pipeline_time_series_to_spectral_connectivity(main_path)

    spectral_workflow.inputs.inputnode.is_sensor_space = False
    
    main_workflow.connect(inv_sol_workflow, 'inv_solution.ts_file',
                          spectral_workflow, 'inputnode.ts_file')

    main_workflow.connect(inv_sol_workflow, 'inv_solution.labels',
                          spectral_workflow, 'inputnode.labels_file')

    main_workflow.connect(infosource, ('freq_band_name', get_freq_band),
                          spectral_workflow, 'inputnode.freq_band')

    main_workflow.connect(preproc_workflow, 'preproc.sfreq',
                          spectral_workflow, 'inputnode.sfreq')

    spectral_workflow.inputs.inputnode.epoch_window_length = 3.0

    graph_den_pipe = create_pipeline_conmat_to_graph_density(main_path,con_den = 0.1,mod = mod,plot = True)

    main_workflow.connect(spectral_workflow, 'spectral.conmat_file',
                          graph_den_pipe, 'compute_net_List.Z_cor_mat_file')
    if mod:                       
        graph_den_pipe.inputs.community_rada.optim_seq = radatools_optim 
        main_workflow.connect(inv_sol_workflow, 'inv_solution.label_names',
                          graph_den_pipe, 'plot_igraph_modules_rada.labels_file')
        main_workflow.connect(inv_sol_workflow, 'inv_solution.label_coords',
                          graph_den_pipe, 'plot_igraph_modules_rada.coords_file')
                          
    # run pipeline
    main_workflow.write_graph(graph2use='colored')  # colored
    main_workflow.config['execution'] = {'remove_unnecessary_outputs': 'false'}
    
    main_workflow.run()
    #main_workflow.run(plugin='MultiProc', plugin_args={'n_procs': 8})
