# -*- coding: utf-8 -*-
"""
Description: 

Wraps spectral connectivity function of MNE, as well as plot_circular_connectivity

"""
import nipype.pipeline.engine as pe
from nipype.interfaces.utility import IdentityInterface
#,Function

from neuropype_ephy.interfaces.mne.spectral import  SpectralConn,PlotSpectralConn
from neuropype_ephy.nodes.ts_tools import SplitWindows

### to modify and add in "Nodes"
#from neuropype_ephy.spectral import  filter_adj_plot_mat

def create_pipeline_time_series_to_spectral_connectivity( main_path,sfreq, pipeline_name = "ts_to_conmat",con_method = "coh", multicon = False, export_to_matlab = False, temporal_windows = []):
    
    if isinstance(main_path,string)
		
		pipeline = pe.Workflow(name= pipeline_name)
		pipeline.base_dir = main_path
		
	else:
		
		pipeline = main_path
		
    inputnode = pe.Node(IdentityInterface(fields=['ts_file','freq_band','labels_file','epoch_window_length','is_sensor_space','index']), name='inputnode')
     
    if multicon == True:
        
        print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Multiple trials $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
        #### spectral
        spectral = pe.MapNode(interface = SpectralConn(), name = "spectral",iterfield = ['ts_file','index'])
        
        spectral.inputs.con_method = con_method  
        spectral.inputs.export_to_matlab = export_to_matlab
        
        spectral.inputs.sfreq = sfreq
        
        pipeline.connect(inputnode, 'ts_file', spectral, 'ts_file')
        pipeline.connect(inputnode, 'freq_band', spectral, 'freq_band')
        pipeline.connect(inputnode, 'index', spectral, 'index')
        pipeline.connect(inputnode, 'epoch_window_length', spectral, 'epoch_window_length')

        #### plot spectral
        plot_spectral = pe.MapNode(interface = PlotSpectralConn(), name = "plot_spectral",iterfield = ['conmat_file'])
        
        pipeline.connect(inputnode,  'labels_file',plot_spectral,'labels_file')
        pipeline.connect(inputnode,  'is_sensor_space',plot_spectral,'is_sensor_space')
        
        pipeline.connect(spectral, "conmat_file",    plot_spectral, 'conmat_file')
        
    else:

        if len(temporal_windows) != 0:

            print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Multiple windows $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
            print n_windows
            ### win_ts
            ##### 
            win_ts = pe.Node(interface = SplitWindows(), name = "win_ts")
            win_ts.inputs.n_windows = n_windows
            
            pipeline.connect(inputnode,'ts_file',win_ts,'ts_file')
            
            win_ts = pe.Node(interface = SplitWindows(), name = "win_ts")
            win_ts.inputs.n_windows = n_windows
            
            
            pipeline.connect(inputnode,'ts_file',win_ts,'ts_file')
            
            #### spectral
            spectral = pe.MapNode(interface = SpectralConn(), name = "spectral",iterfield = ['ts_file'])
            spectral.inputs.con_method = con_method  
            spectral.inputs.export_to_matlab = False
            spectral.inputs.sfreq = sfreq 
            #spectral.inputs.freq_band = freq_band 
            main_workflow.connect(win_ts, 'win_ts_files', spectral, 'ts_file')
            main_workflow.connect(infosource, ('freq_band_name',get_freq_band),spectral, 'freq_band')
            
            ##### plot spectral
            plot_spectral = pe.MapNode(interface = PlotSpectralConn(), name = "plot_spectral", iterfield = ['conmat_file'])
            
            main_workflow.connect(datasource,  'channel_names_file',plot_spectral,'labels_file')
            
            #plot_spectral.inputs.is_sensor_space = False
            main_workflow.connect(spectral, "conmat_file",    plot_spectral, 'conmat_file')


            
        else:
                
            #### spectral
            spectral = pe.Node(interface = SpectralConn(), name = "spectral")
            
            spectral.inputs.con_method = con_method  
            spectral.inputs.export_to_matlab = export_to_matlab            
            spectral.inputs.sfreq = sfreq
        
            pipeline.connect(inputnode, 'ts_file', spectral, 'ts_file')
            pipeline.connect(inputnode, 'freq_band', spectral, 'freq_band')
            pipeline.connect(inputnode, 'epoch_window_length', spectral, 'epoch_window_length')

            #### plot spectral
            plot_spectral = pe.Node(interface = PlotSpectralConn(), name = "plot_spectral")
            
            pipeline.connect(inputnode,  'labels_file',plot_spectral,'labels_file')
            pipeline.connect(inputnode,  'is_sensor_space',plot_spectral,'is_sensor_space')
            
            pipeline.connect(spectral, "conmat_file",    plot_spectral, 'conmat_file')
            
            

    return pipeline
    