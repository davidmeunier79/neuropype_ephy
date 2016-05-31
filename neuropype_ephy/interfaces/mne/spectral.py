# -*- coding: utf-8 -*-

"""
Definition of nodes for computing and plotting spectral connectivity
"""
import numpy as np
import os
import pickle

from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec, isdefined
    
from nipype.utils.filemanip import split_filename as split_f
    
############################################################################################### SpectralConn #####################################################################################################

from neuropype_ephy.spectral import compute_and_save_spectral_connectivity

class SpectralConnInputSpec(BaseInterfaceInputSpec):
    
    ts_file = traits.File(exists=True, desc='nodes * time series in .npy format', mandatory=True)
    
    sfreq = traits.Float(desc='sampling frequency', mandatory=True)
    
    freq_band = traits.List(traits.Float(exists=True), desc='frequency bands', mandatory=True)
    
    con_method = traits.Enum("coh","imcoh","plv","pli","wpli","pli2_unbiased","ppc","cohy","wpli2_debiased",desc='metric computed on time series for connectivity')
    
    epoch_window_length = traits.Float(desc='epoched data', mandatory=False)
    
    export_to_matlab = traits.Bool(False, desc='If conmat is exported to .mat format as well',usedefault = True)
    
class SpectralConnOutputSpec(TraitedSpec):
    
    conmat_file = File(exists=True, desc="spectral connectivty matrix in .npy format")
    
class SpectralConn(BaseInterface):
    
    """
    Description:
    
    Compute spectral connectivity in a given frequency bands
    
    Inputs:
    
    ts_file 
        type = File , exists=True, desc='nodes * time series in .npy format', mandatory=True
    
    sfreq 
        type = Float, desc='sampling frequency', mandatory=True
    
    freq_band 
        type = List(Float) , exists=True, desc='frequency bands', mandatory=True
    
    con_method 
        type = Enum("coh","imcoh","plv","pli","wpli","pli2_unbiased","ppc","cohy","wpli2_debiased") , desc='metric computed on time series for connectivity'
        
    epoch_window_length 
        type = Float, desc='epoched data', mandatory=False
    
    export_to_matlab 
        type = Bool, default = False, desc='If conmat is exported to .mat format as well',usedefault = True
   
    Outputs:
    
    conmat_file 
        type = File, exists=True, desc="spectral connectivty matrix in .npy format"
    
    
    """
    input_spec = SpectralConnInputSpec
    output_spec = SpectralConnOutputSpec

    def _run_interface(self, runtime):
                
        print 'in SpectralConn'
        
        ts_file = self.inputs.ts_file
        sfreq = self.inputs.sfreq
        freq_band = self.inputs.freq_band
        con_method = self.inputs.con_method
        epoch_window_length = self.inputs.epoch_window_length
        export_to_matlab = self.inputs.export_to_matlab
        
        if epoch_window_length == traits.Undefined:
            data = np.load(ts_file)
        else:
            raw_data = np.load(ts_file)
            nb_splits = raw_data.shape[1] // (epoch_window_length * sfreq)
            reste = raw_data.shape[1] % int(epoch_window_length * sfreq)
            if reste != 0:
                raw_data = raw_data[:,:-reste]
            print "epoching data with {}s by window, resulting in {} epochs (rest = {})".format(epoch_window_length,nb_splits,reste)
            data = np.array(np.array_split(raw_data,nb_splits,axis = 1))
        
        self.conmat_file = compute_and_save_spectral_connectivity(data = data,con_method = con_method,sfreq=sfreq, fmin= freq_band[0], fmax=freq_band[1],export_to_matlab = export_to_matlab)
        
        return runtime
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["conmat_file"] = self.conmat_file
        
        return outputs
        
############################################################################################### PlotSpectralConn #####################################################################################################

from neuropype_ephy.spectral import plot_circular_connectivity

class PlotSpectralConnInputSpec(BaseInterfaceInputSpec):
    
    conmat_file = traits.File(exists=True, desc='connectivity matrix in .npy format', mandatory=True)
    
    is_sensor_space = traits.Bool(True, desc = 'if True uses labels as returned from mne', usedefault = True)
    
    vmin = traits.Float(0.3, desc='min scale value', usedefault = True)
    
    vmax = traits.Float(1.0, desc='max scale value', usedefault = True)
    
    nb_lines = traits.Int(200, desc='nb lines kept in the representation', usedefault = True)
    
    labels_file = traits.File(desc='list of labels associated with nodes')
    
class PlotSpectralConnOutputSpec(TraitedSpec):
    
    plot_conmat_file = File(exists=True, desc="plot spectral connectivity matrix in .eps format")
    
class PlotSpectralConn(BaseInterface):
    
    """
    Description:
    
    Plot connectivity matrix using mne plot_circular_connectivity function
    
    Inputs:
    
    conmat_file 
        type = File, exists=True, desc='connectivity matrix in .npy format', mandatory=True
    
    is_sensor_space 
        type = Bool, default = True, desc = 'if True uses labels as returned from mne', usedefault = True
    
    vmin 
        type = Float, default = 0.3, desc='min scale value', usedefault = True
    
    vmax
        type = Float, default = 1.0, desc='max scale value', usedefault = True
    
    nb_lines 
        type = Int, default = 200, desc='nb lines kept in the representation', usedefault = True
    
    labels_file 
        type = File, desc='list of labels associated with nodes'
    
    Outputs:
    
    plot_conmat_file
        type = File, exists=True, desc="plot spectral connectivity matrix in .eps format"
        
    
    """
    
    input_spec = PlotSpectralConnInputSpec
    output_spec = PlotSpectralConnOutputSpec

    def _run_interface(self, runtime):
                
        print 'in PlotSpectralConn'
        
        conmat_file = self.inputs.conmat_file
        vmin = self.inputs.vmin
        vmax = self.inputs.vmax
        nb_lines = self.inputs.nb_lines
        is_sensor_space = self.inputs.is_sensor_space
        labels_file =self.inputs.labels_file
        
        
        ### reading matrix and base filename from conmat_file
        path,fname,ext = split_f(conmat_file)    
        print fname
        
        conmat = np.load(conmat_file)
        print conmat.shape
        
        assert conmat.ndim == 2, "Warning, conmat should be 2D matrix , ndim = {}".format(conmat.ndim)
        assert conmat.shape[0] == conmat.shape[1], "Warning, conmat should be a squared matrix , {} != {}".format(conmat.shape[0],conmat.shape[1])
        
            
        if isdefined(labels_file):
            
            if is_sensor_space:
                label_names = [line.strip() for line in open(labels_file)]
                node_order  = label_names
                node_colors = None
            
            else:
                labels = []
                with open(labels_file, "rb") as f:
                    for _ in range(pickle.load(f)):
                        labels.append(pickle.load(f))
                
                print labels
#                0/0
                # read colors
                node_colors = [label.color for label in labels]    
                
                # reorder the labels based on their location in the left hemi
                label_names = [label.name for label in labels]
                lh_labels = [name for name in label_names if name.endswith('lh')]

                # Get the y-location of the label
                label_ypos = list()
                for name in lh_labels:
                    idx = label_names.index(name)
                    ypos = np.mean(labels[idx].pos[:, 1])
                    label_ypos.append(ypos)

                try:
                    idx = label_names.index('Brain-Stem')
                    ypos = np.mean(labels[idx].pos[:, 1])
                    lh_labels.append('Brain-Stem')
                    label_ypos.append(ypos)
                except ValueError:
                    pass
                        
                # Reorder the labels based on their location
                lh_labels = [label for (yp, label) in sorted(zip(label_ypos, lh_labels))]
                
                # For the right hemi
                rh_labels = [label[:-2] + 'rh' for label in lh_labels if label != 'Brain-Stem']

                # Save the plot order 
                node_order = list()
                node_order.extend(lh_labels[::-1])  # reverse the order
                node_order.extend(rh_labels)             

        else:
            label_names = range(conmat.shape[0])
            node_order  = label_names
            node_colors = None
           
        print label_names
        print node_order
#        0/0
        self.plot_conmat_file = plot_circular_connectivity(conmat,label_names,node_colors,node_order, vmin,vmax ,nb_lines, fname)

        return runtime
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["plot_conmat_file"] = self.plot_conmat_file
        
        return outputs
        
