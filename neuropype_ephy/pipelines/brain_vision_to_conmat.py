# -*- coding: utf-8 -*-

import nipype.pipeline.engine as pe

from nipype.interfaces.utility import Function
from nipype.interfaces.utility import IdentityInterface


from neuropype_ephy.interfaces.mne.spectral import  SpectralConn,PlotSpectralConn

#from neuropype_ephy.spectral import  multiple_spectral_proc

from neuropype_ephy.nodes.import_data import ImportBrainVisionAscii

from neuropype_ephy.nodes.ts_tools import SplitWindows

#from neuropype_ephy.spectral import split_win_ts

###TODO
#from neuropype_ephy.nodes.? import filter_adj_plot_mat

def create_pipeline_brain_vision_ascii_to_spectral_connectivity(main_path,pipeline_name="brain_vision_to_conmat", con_method = "coh", sample_size = 512, sep_label_name = "", sfreq = 512,filter_spectral = True, k_neigh = 3, n_windows = []):
    
    """
    Description:
    
    Create pipeline from intraEEG times series in ascii format exported out of BrainVision, split txt and compute spectral connectivity.
    Possibly also filter out connections between "adjacent" contacts (on the same electrode)
    
    
    """
    pipeline = pe.Workflow(name=pipeline_name )
    pipeline.base_dir = main_path
    
    
    inputnode = pe.Node(interface = IdentityInterface(fields=['txt_file','freq_band']), name='inputnode')
    
    
    #### convert
    #split_ascii = pe.Node(interface = Function(input_names = ["sample_size","txt_file","sep_label_name"],output_names = ["splitted_ts_file","elec_names_file"],function = split_txt),name = 'split_ascii')
    
    split_ascii = pe.Node(interface = ImportBrainVisionAscii(),name = 'split_ascii')
    
    split_ascii.inputs.sample_size = sample_size
    split_ascii.inputs.sep_label_name = sep_label_name
    
    pipeline.connect(inputnode, 'txt_file',split_ascii,'txt_file')

    if len(n_windows) == 0:
            
        #### spectral
        
        spectral = pe.Node(interface = SpectralConn(), name = "spectral")
        
        #spectral = pe.Node(interface = Function(input_names = ["ts_file","sfreq","freq_band","con_method"],
        #                                        output_names = "conmat_file",
        #                                        function = spectral_proc),name = "spectral")
        
        spectral.inputs.con_method = con_method    
        spectral.inputs.sfreq = sfreq
        
        pipeline.connect(inputnode, 'freq_band', spectral, 'freq_band')
        
        pipeline.connect(split_ascii, 'splitted_ts_file', spectral, 'ts_file')

        #### plot spectral
        plot_spectral = pe.Node(interface = PlotSpectralConn(), name = "plot_spectral")
        
        #plot_spectral = pe.Node(interface = Function(input_names = ["conmat_file","labels_file","nb_lines","vmin","vmax"],
                                                    #output_names = "plot_conmat_file",
                                                    #function = plot_circular_connectivity), name = "plot_spectral")
        
        # plot_spectral.inputs.labels_file = MEG_elec_names_file AP 021015
        plot_spectral.inputs.nb_lines = 200
        plot_spectral.inputs.vmin = 0.3
        plot_spectral.inputs.vmax = 1.0
        
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
            plot_filter_spectral = pe.Node(interface = Function(input_names = ["conmat_file","labels_file","nb_lines","vmin","vmax"],
                                                        output_names = "plot_conmat_file",
                                                        function = plot_circular_connectivity), name = "plot_filter_spectral_" + str(k_neigh))
            
            # plot_spectral.inputs.labels_file = MEG_elec_names_file AP 021015
            plot_filter_spectral.inputs.nb_lines = 50
            
            plot_filter_spectral.inputs.vmin = 0.3
            plot_filter_spectral.inputs.vmax = 1.0
        
        
            pipeline.connect(split_ascii,  'elec_names_file',plot_filter_spectral,'labels_file')
            pipeline.connect(filter_spectral, "filtered_conmat_file",    plot_filter_spectral, 'conmat_file')
            
    else:
            
        ### win_ts
        ##### 
        win_ts = pe.Node(interface = SplitWindows(), name = "win_ts")
        
        win_ts.inputs.n_windows = n_windows
        
        pipeline.connect(split_ascii,'splitted_ts_file',win_ts,'ts_file')
                 
        #win_ts = pe.Node(interface = Function(input_names = ["splitted_ts_file","n_windows"],output_names = ["win_splitted_ts_files"],function = split_win_ts), name = "win_ts")
        
        #win_ts.inputs.n_windows = n_windows
        
        #pipeline.connect(split_ascii,'splitted_ts_file',win_ts,'splitted_ts_file')
                 
                 
        spectral = pe.MapNode(interface = SpectralConn(), iterfield = ['ts_file'], name = "spectral")
        
        
        spectral.inputs.con_method = con_method    
        spectral.inputs.sfreq = sfreq
        
        #spectral.inputs.epoch_window_length = epoch_window_length
        pipeline.connect(win_ts, 'win_ts_files', spectral, 'ts_file')
        pipeline.connect(inputnode,'freq_band', spectral, 'freq_band')
        
        ##### spectral
        #spectral = pe.Node(interface = Function(input_names = ["ts_file","sfreq","freq_band","con_method"],
                                                #output_names = "conmat_files",
                                                #function = multiple_spectral_proc),name = "spectral")
        
        #spectral.inputs.con_method = con_method    
        #spectral.inputs.sfreq = sfreq
        
        #pipeline.connect(split_ascii, 'splitted_ts_file', spectral, 'ts_file')

    return pipeline
    
    
#def create_pipeline_brain_vision_ascii_to_multiwin_spectral_connectivity(main_path,t_windows, con_method = "coh", sample_size = 512, sep_label_name = "", sfreq = 512):

    #pipeline = pe.Workflow(name="brain_vision_to_multiwin_conmat")
    #pipeline.base_dir = main_path
    
    ##### convert
    #split_ascii_multiwin = pe.Node(interface = Function(input_names = ["sample_size","txt_file","sep_label_name","t_windows"],output_names = ["splitted_ts_file","elec_names_file"],function = split_txt_multiwin),name = 'split_ascii_multiwin')
    
    #split_ascii_multiwin.inputs.t_windows = t_windows
    #split_ascii_multiwin.inputs.sample_size = sample_size
    #split_ascii_multiwin.inputs.sep_label_name = sep_label_name
    
    ###### spectral
    ##spectral = pe.Node(interface = Function(input_names = ["ts_file","sfreq","freq_band","con_method"],
                                            ##output_names = "conmat_file",
                                            ##function = spectral_proc),name = "spectral")
    
    ##spectral.inputs.con_method = con_method    
    ##spectral.inputs.sfreq = sfreq
    
    ##pipeline.connect(split_ascii, 'splitted_ts_file', spectral, 'ts_file')

    ###### plot spectral
    ##plot_spectral = pe.Node(interface = Function(input_names = ["conmat_file","labels_file"],
                                                 ##output_names = "plot_conmat_file",
                                                 ##function = plot_circular_connectivity), name = "plot_spectral")
    
    ### plot_spectral.inputs.labels_file = MEG_elec_names_file AP 021015
    ##plot_spectral.inputs.nb_lines = 200
    
    ##pipeline.connect(split_ascii,  'elec_names_file',plot_spectral,'labels_file')
    ##pipeline.connect(spectral, "conmat_file",    plot_spectral, 'conmat_file')
    
    #return pipeline
    
    
        
    
    
    

    

