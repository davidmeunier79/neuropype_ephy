# -*- coding: utf-8 -*-

"""
Definition of nodes for computing reordering and plotting coclass_matrices
"""
import numpy as np
import os

from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec, isdefined
    
############################################################################################### SpectralConn #####################################################################################################

from neuropype_ephy.spectral import spectral_proc,epoched_spectral_proc

class SpectralConnInputSpec(BaseInterfaceInputSpec):
    
    ts_file = traits.File(exists=True, desc='nodes * time series in .npy format', mandatory=True)
    
    sfreq = traits.Float(desc='sampling frequency', mandatory=True)
    
    freq_band = traits.List(traits.Float(exists=True), desc='frequency bands', mandatory=True)
    
    con_method = traits.Enum("coh","imcoh","plv","pli",desc='metric computed on time series for connectivity')
    
    epoch_window_length = traits.Float(desc='epoched data', mandatory=False)
    
class SpectralConnOutputSpec(TraitedSpec):
    
    conmat_file = File(exists=True, desc="spectral connectivty matrix in .npy format")
    
class SpectralConn(BaseInterface):
    
    """
    Compute spectral connectivity in a given frequency bands
    """
    input_spec = SpectralConnInputSpec
    output_spec = SpectralConnOutputSpec

    def _run_interface(self, runtime):
                
        print 'in prepare_coclass'
        ts_file = self.inputs.ts_file
        sfreq = self.inputs.sfreq
        freq_band = self.inputs.freq_band
        con_method = self.inputs.con_method
        epoch_window_length = self.inputs.epoch_window_length
        
        if epoch_window_length == traits.Undefined:
            self.conmat_file = spectral_proc(ts_file,sfreq,freq_band,con_method)
        
        else:
            self.conmat_file = epoched_spectral_proc(ts_file,sfreq,freq_band,con_method,epoch_window_length)

        return runtime
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["conmat_file"] = self.conmat_file
        
        return outputs
        
