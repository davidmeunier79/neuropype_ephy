# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 17:06:36 2016

@author: pasca
"""

# -*- coding: utf-8 -*-


import matplotlib
matplotlib.use('PS')


def get_ext_file(raw_file):
    from nipype.utils.filemanip import split_filename as split_f

    subj_path, basename, ext = split_f(raw_file)

    print raw_file
    is_ds = False
    if ext is 'ds':
        is_ds = True
        return is_ds
    elif ext is 'fif':
        return is_ds
    else:
        raise RuntimeError('only fif and ds file format!!!')


# is_ICA=True                 => apply ICA to automatically remove ECG and EoG
#                                artifacts
# is_set_ICA_components=False => specify all subject_ids and sessions
# is_set_ICA_components=True  => specify the dataset for we want to recompute
#                               the ICA
# in Elekta data, ICA routine automatically looks for EEG61, EEG62
def create_pipeline_preproc_meeg(main_path,
                                 pipeline_name='preproc_meeg',
                                 data_type='fif',
                                 l_freq=1, h_freq=150, down_sfreq=300,
                                 is_ICA=True, variance=0.95,
                                 ECG_ch_name='', EoG_ch_name='',
                                 is_set_ICA_components=False,
                                 n_comp_exclude=[],
                                 is_sensor_space=True):

    from neuropype_ephy.preproc import preprocess_fif_to_ts
    from neuropype_ephy.preproc import preprocess_ICA_fif_to_ts
    from neuropype_ephy.preproc import preprocess_set_ICA_comp_fif_to_ts
    from nipype.interfaces.utility import IdentityInterface, Function
    from neuropype_ephy.import_ctf import convert_ds_to_raw_fif

    import nipype
    print nipype.__version__

    import nipype.pipeline.engine as pe

    pipeline = pe.Workflow(name=pipeline_name)
    pipeline.base_dir = main_path

    print '*** main_path -> %s' % main_path + ' ***'
    print '*** is_sensor_space -> %s ***' % is_sensor_space

    # define the inputs of the pipeline
    inputnode = pe.Node(IdentityInterface(fields=['raw_file']),
                        name='inputnode')

    if data_type is 'ds':
        convert = pe.Node(interface=Function(input_names=['ds_file'],
                                             output_names=['raw_fif_file'],
                                             function=convert_ds_to_raw_fif),
                          name='convert_ds')

        pipeline.connect(inputnode, 'raw_file', convert, 'ds_file')

    # preprocess
    if is_ICA:
        if is_set_ICA_components:
            preproc = pe.Node(interface=Function(input_names=['fif_file',
                                                              'n_comp_exclude',
                                                              'l_freq',
                                                              'h_freq',
                                                              'down_sfreq',
                                                              'is_sensor_space'],
                                                 output_names=['out_file',
                                                               'channel_coords_file',
                                                               'channel_names_file',
                                                               'sfreq'],
                                                 function=preprocess_set_ICA_comp_fif_to_ts),
                              name='preproc')
            preproc.inputs.n_comp_exclude = n_comp_exclude
        else:
            preproc = pe.Node(interface=Function(input_names=['fif_file',
                                                              'ECG_ch_name',
                                                              'EoG_ch_name',
                                                              'l_freq',
                                                              'h_freq',
                                                              'down_sfreq',
                                                              'variance',
                                                              'is_sensor_space',
                                                              'data_type'],
                                                 output_names=['out_file',
                                                               'channel_coords_file',
                                                               'channel_names_file',
                                                               'sfreq'],
                                                 function=preprocess_ICA_fif_to_ts),
                              name='preproc')
            preproc.inputs.ECG_ch_name = ECG_ch_name
            preproc.inputs.EoG_ch_name = EoG_ch_name
            preproc.inputs.data_type = data_type
            preproc.inputs.variance = variance

    else:
        preproc = pe.Node(interface=Function(input_names=['fif_file',
                                                          'l_freq',
                                                          'h_freq',
                                                          'down_sfreq'],
                                             output_names=['out_file',
                                                           'channel_coords_file',
                                                           'channel_names_file',
                                                           'sfreq'],
                                             function=preprocess_fif_to_ts),
                          name='preproc')

    preproc.inputs.is_sensor_space = is_sensor_space
    preproc.inputs.l_freq = l_freq
    preproc.inputs.h_freq = h_freq
    preproc.inputs.down_sfreq = down_sfreq

    if data_type is 'ds':
        pipeline.connect(convert, 'raw_fif_file', preproc, 'fif_file')
    elif data_type is 'fif':
        pipeline.connect(inputnode, 'raw_file', preproc, 'fif_file')

    return pipeline

# TODO remove this TEST code
if __name__ == '__main__':
    main_path = '/home/karim/Documents/Fanny'
    data_path = main_path
    subject_ids = ['S01']
    sessions = ['repos_1']

    # create main workflow
    main_workflow = pe.Workflow(name='wf_preproc_meeg')
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

    '''
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                   outfields=['fif_file']),
                         name='datasource')
    datasource.inputs.base_directory = data_path
    datasource.inputs.template = '%s/%s*%s'
    datasource.inputs.template_args = dict(fif_file=[['subject_id',
                                                      'subject_id',
                                                      "_tsss_mc.fif"]])
    datasource.inputs.sort_filelist = True
    '''
    main_workflow.connect(infosource, 'subject_id', datasource, 'subject_id')
    main_workflow.connect(infosource, 'sess_index', datasource, 'sess_index')

#    preproc_workflow = create_pipeline_preproc_meeg(main_path,
#                                                    is_set_ICA_components=True,
#                                                    n_comp_exclude={'S01': [[0,5,53], [1,19,28]], 
#                                                                    'S02': [[1,10,18,78],[6,33,52]]})

    preproc_workflow = create_pipeline_preproc_meeg(main_path, data_type='ds')

    main_workflow.connect(datasource, 'raw_file',
                          preproc_workflow, 'inputnode.raw_file')

    # run pipeline
    main_workflow.write_graph(graph2use='colored')  # colored
    main_workflow.config['execution'] = {'remove_unnecessary_outputs': 'false'}
    main_workflow.run(plugin='MultiProc', plugin_args={'n_procs': 8})
