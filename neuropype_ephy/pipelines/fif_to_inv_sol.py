# -*- coding: utf-8 -*-

import os
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio

from nipype.interfaces.utility import IdentityInterface
from neuropype_ephy.interfaces.mne.LF_computation import LFComputation
from neuropype_ephy.interfaces.mne.Inverse_solution import NoiseCovariance
from neuropype_ephy.pipelines.preproc_meeg import create_pipeline_preproc_meeg


def create_pipeline_source_reconstruction(main_path, sbj_dir,
                                          pipeline_name='inv_sol_pipeline',
                                          spacing='ico-5',
                                          aseg=False,
                                          aseg_labels=[],
                                          noise_cov_fname=None):

    pipeline = pe.Workflow(name=pipeline_name)
    pipeline.base_dir = main_path

    inputnode = pe.Node(IdentityInterface(fields=['sbj_id', 'raw']),
                        name='inputnode')

    LF_computation = pe.Node(interface=LFComputation(), name='LF_computation')

    LF_computation.inputs.sbj_dir = sbj_dir
    LF_computation.inputs.spacing = spacing
    if aseg:
        LF_computation.inputs.aseg = aseg
        LF_computation.inputs.aseg_labels = aseg_labels

    pipeline.connect(inputnode, 'sbj_id', LF_computation, 'sbj_id')

    pipeline.connect(inputnode, ('raw', get_raw_info),
                     LF_computation, 'raw_info')

    create_noise_cov = pe.Node(interface=NoiseCovariance(),
                               name="create_noise_cov")

    if noise_cov_fname is not None:
        create_noise_cov.inputs.cov_fname_in = noise_cov_fname

    pipeline.connect(inputnode, 'raw', create_noise_cov, 'raw')

    return pipeline


def get_raw_info(raw):
    return raw.info

if __name__ == '__main__':
    main_path = '/home/karim/Documents/Fanny'
    data_path = main_path
    sbj_dir = os.path.join(main_path, 'FSF')
    subject_ids = ['S01']
    sessions = ['repos_1']  # 'repos_2'

    noise_cov_fname = os.path.join(main_path, 'Big_Noise-cov.fif')

    aseg_labels = ['Left-Accumbens-area',
                   'Left-Amygdala',
                   'Left-Caudate',
                   'Left-Hippocampus',
                   'Left-Pallidum',
                   'Left-Putamen',
                   'Left-Thalamus-Proper',
                   'Left-Cerebellum-Cortex',
                   'Brain-Stem',
                   'Right-Accumbens-area',
                   'Right-Amygdala',
                   'Right-Caudate',
                   'Right-Hippocampus',
                   'Right-Pallidum',
                   'Right-Putamen',
                   'Right-Thalamus-Proper',
                   'Right-Cerebellum-Cortex']

    # create main workflow
    main_workflow = pe.Workflow(name='wf_LF_computation')
    main_workflow.base_dir = main_path

    # info source
    infosource = pe.Node(interface=IdentityInterface(fields=['subject_id',
                                                             'sess_index']),
                         name="infosource")
    infosource.iterables = [('subject_id', subject_ids),
                            ('sess_index', sessions)]

    # data source
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id',
                                                             'sess_index'],
                                                   outfields=['raw_file']),
                         name='datasource')
    datasource.inputs.base_directory = data_path
    datasource.inputs.template = '%s/%s%s'
    datasource.inputs.template_args = dict(raw_file=[['subject_id',
                                                     'sess_index', ".ds"]])
    datasource.inputs.sort_filelist = True

    main_workflow.connect(infosource, 'subject_id', datasource, 'subject_id')
    main_workflow.connect(infosource, 'sess_index', datasource, 'sess_index')

    preproc_workflow = create_pipeline_preproc_meeg(main_path,
                                                    is_sensor_space=False,
                                                    data_type='ds')

    main_workflow.connect(datasource, 'raw_file',
                          preproc_workflow, 'inputnode.raw_file')

    inv_sol_workflow = create_pipeline_source_reconstruction(main_path,
                                                             sbj_dir,
                                                             noise_cov_fname=noise_cov_fname)
#    inv_sol_workflow = create_pipeline_source_reconstruction(main_path,
#                                                             sbj_dir)

    main_workflow.connect(infosource, 'subject_id',
                          inv_sol_workflow, 'inputnode.sbj_id')

    main_workflow.connect(preproc_workflow, 'preproc.out_file',
                          inv_sol_workflow, 'inputnode.raw')

    # run pipeline
    main_workflow.write_graph(graph2use='colored')  # colored
    main_workflow.config['execution'] = {'remove_unnecessary_outputs': 'false'}
    main_workflow.run(plugin='MultiProc', plugin_args={'n_procs': 8})
