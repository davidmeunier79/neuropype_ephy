# -*- coding: utf-8 -*-

import nipype.pipeline.engine as pe
from nipype.interfaces.utility import Function

from neuropype_ephy.import_txt import split_txt,split_txt_multiwin

from neuropype_ephy.spectral import  epoched_spectral_proc,spectral_proc
from neuropype_ephy.spectral import  plot_circular_connectivity,filter_adj_plot_mat

def create_pipeline_brain_vision_ascii_to_spectral_connectivity(main_path,con_method = "coh", sample_size = 512, sep_label_name = "", sfreq = 512,filter_spectral = True, k_neigh = 3):

    pipeline = pe.Workflow(name="brain_vision_to_conmat")
    pipeline.base_dir = main_path
    
    #### convert
    split_ascii = pe.Node(interface = Function(input_names = ["sample_size","txt_file","sep_label_name"],output_names = ["splitted_ts_file","elec_names_file"],function = split_txt),name = 'split_ascii')
    split_ascii.inputs.sample_size = sample_size
    split_ascii.inputs.sep_label_name = sep_label_name
    
    #### spectral
    spectral = pe.Node(interface = Function(input_names = ["ts_file","sfreq","freq_band","con_method"],
                                            output_names = "conmat_file",
                                            function = spectral_proc),name = "spectral")
    
    spectral.inputs.con_method = con_method    
    spectral.inputs.sfreq = sfreq
    
    pipeline.connect(split_ascii, 'splitted_ts_file', spectral, 'ts_file')

    #### plot spectral
    plot_spectral = pe.Node(interface = Function(input_names = ["conmat_file","labels_file"],
                                                 output_names = "plot_conmat_file",
                                                 function = plot_circular_connectivity), name = "plot_spectral")
    
    # plot_spectral.inputs.labels_file = MEG_elec_names_file AP 021015
    plot_spectral.inputs.nb_lines = 200
    
    pipeline.connect(split_ascii,  'elec_names_file',plot_spectral,'labels_file')
    pipeline.connect(spectral, "conmat_file",    plot_spectral, 'conmat_file')
    
    
    if filter_spectral == True:
            
        ### filter spectral
        filter_spectral = pe.Node(interface = Function(input_names = ["conmat_file","labels_file","sep_label_name","k_neigh"], 
                                                    output_names = "filtered_conmat_file", 
                                                    function = filter_adj_plot_mat), name = "filter_spectral_" + str(k_neigh))
        filter_spectral.inputs.sep_label_name = sep_label_name
        filter_spectral.inputs.k_neigh = k_neigh
        
        pipeline.connect(split_ascii,  'elec_names_file',filter_spectral,'labels_file')
        pipeline.connect(spectral, "conmat_file",    filter_spectral, 'conmat_file')
        
        
        
        #### plot filter_spectral
        plot_filter_spectral = pe.Node(interface = Function(input_names = ["conmat_file","labels_file"],
                                                    output_names = "plot_conmat_file",
                                                    function = plot_circular_connectivity), name = "plot_filter_spectral_" + str(k_neigh))
        
        # plot_spectral.inputs.labels_file = MEG_elec_names_file AP 021015
        plot_filter_spectral.inputs.nb_lines = 50
        
        pipeline.connect(split_ascii,  'elec_names_file',plot_filter_spectral,'labels_file')
        pipeline.connect(filter_spectral, "filtered_conmat_file",    plot_filter_spectral, 'conmat_file')
        
    
    return pipeline
    
    
def create_pipeline_brain_vision_ascii_to_multiwin_spectral_connectivity(main_path,t_windows, con_method = "coh", sample_size = 512, sep_label_name = "", sfreq = 512):

    pipeline = pe.Workflow(name="brain_vision_to_multiwin_conmat")
    pipeline.base_dir = main_path
    
    #### convert
    split_ascii_multiwin = pe.Node(interface = Function(input_names = ["sample_size","txt_file","sep_label_name","t_windows"],output_names = ["splitted_ts_file","elec_names_file"],function = split_txt_multiwin),name = 'split_ascii_multiwin')
    
    split_ascii_multiwin.inputs.t_windows = t_windows
    split_ascii_multiwin.inputs.sample_size = sample_size
    split_ascii_multiwin.inputs.sep_label_name = sep_label_name
    
    ##### spectral
    #spectral = pe.Node(interface = Function(input_names = ["ts_file","sfreq","freq_band","con_method"],
                                            #output_names = "conmat_file",
                                            #function = spectral_proc),name = "spectral")
    
    #spectral.inputs.con_method = con_method    
    #spectral.inputs.sfreq = sfreq
    
    #pipeline.connect(split_ascii, 'splitted_ts_file', spectral, 'ts_file')

    ##### plot spectral
    #plot_spectral = pe.Node(interface = Function(input_names = ["conmat_file","labels_file"],
                                                 #output_names = "plot_conmat_file",
                                                 #function = plot_circular_connectivity), name = "plot_spectral")
    
    ## plot_spectral.inputs.labels_file = MEG_elec_names_file AP 021015
    #plot_spectral.inputs.nb_lines = 200
    
    #pipeline.connect(split_ascii,  'elec_names_file',plot_spectral,'labels_file')
    #pipeline.connect(spectral, "conmat_file",    plot_spectral, 'conmat_file')
    
    return pipeline
    
    
        
    
    
    

    

