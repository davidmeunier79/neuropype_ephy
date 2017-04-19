# Created on Tue Apr 26 17:06:36 2016
# @author: pasca


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
                                 reject=None,
                                 is_set_ICA_components=False,
                                 n_comp_exclude=[]):

    """
    Description:
    
        Preprocessing pipeline

    Inputs:

        main_path : str
            the main path of the pipeline
        pipeline_name: str (default 'preproc_meeg')
            name of the pipeline
        data_type: str (default 'fif')
            data type: 'fif' or 'ds'
        l_freq: float (default 1)
            low cut-off frequency in Hz
        h_freq: float (default 150)
            high cut-off frequency in Hz
        down_sfreq: float (default 300)
            sampling frequency at which the data are downsampled
        is_ICA : boolean (default True)
            if True apply ICA to automatically remove ECG and EoG artifacts
        variance: float (default 0.95)
            the cumulative percentage of explained variance
        ECG_ch_name: str
            the name of ECG channels
        EoG_ch_name:
            the name of EoG channels
        reject: dict | None
            rejection parameters based on peak-to-peak amplitude. Valid keys
            are 'grad' | 'mag' | 'eeg' | 'eog' | 'ecg'. If reject is None then
            no rejection is done
        is_set_ICA_components: boolean (default False)
            set to True if we had the ICA of the raw data, checked the Report
            and want to exclude some ICA components based on their topographies
            and time series
            if True, we have to fill the dictionary variable n_comp_exclude
        n_comp_exclude: dict
            if is_set_ICA_components=True, it has to be a dict containing for
            each subject and for each session the components to be excluded
            
    Outouts:

        pipeline : instance of Workflow

    """
    
    from neuropype_ephy.interfaces.mne.preproc import PreprocFif
    from neuropype_ephy.interfaces.mne.preproc import CompIca
    from neuropype_ephy.nodes.import_data import ConvertDs2Fif
    from neuropype_ephy.preproc import preprocess_set_ICA_comp_fif_to_ts
    from nipype.interfaces.utility import IdentityInterface, Function

    import nipype
    print nipype.__version__

    import nipype.pipeline.engine as pe

    pipeline = pe.Workflow(name=pipeline_name)
    pipeline.base_dir = main_path

    print '*** main_path -> %s' % main_path + ' ***'

    # define the inputs of the pipeline
    inputnode = pe.Node(IdentityInterface(fields=['raw_file', 'subject_id']),
                        name='inputnode')

    if data_type is 'ds':
        
        ds2fif_node = pe.Node(interface=ConvertDs2Fif(), name='ds2fif')
        pipeline.connect(inputnode, 'raw_file', ds2fif_node, 'ds_file')        

    # preprocess
    preproc_node = pe.Node(interface=PreprocFif(), name='preproc')

    preproc_node.inputs.l_freq = l_freq
    preproc_node.inputs.h_freq = h_freq
    preproc_node.inputs.down_sfreq = down_sfreq
    
    
    if data_type is 'ds':
        pipeline.connect(ds2fif_node, 'fif_file', preproc_node, 'fif_file')
    elif data_type is 'fif':
        pipeline.connect(inputnode, 'raw_file', preproc_node, 'fif_file')
    

    if is_ICA:
        if is_set_ICA_components:
            preproc = pe.Node(interface=Function(input_names=['fif_file',
                                                              'subject_id',
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
            
            pipeline.connect(inputnode, 'subject_id', preproc, 'subject_id')
            
        else:
            
            ica_node = pe.Node(interface=CompIca(), name='ica')
            ica_node.inputs.n_components = variance
            ica_node.inputs.ecg_ch_name = ECG_ch_name
            ica_node.inputs.eog_ch_name = EoG_ch_name

            pipeline.connect(preproc_node, 'fif_file', ica_node, 'fif_file')
            

    return pipeline

