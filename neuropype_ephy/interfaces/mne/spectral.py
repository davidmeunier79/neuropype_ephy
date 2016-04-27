# -*- coding: utf-8 -*-

"""
Definition of nodes for computing reordering and plotting coclass_matrices
"""
import numpy as np
import os

from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec, isdefined
    
############################################################################################### SpectralConn #####################################################################################################


Function(input_names = ["ts_file","sfreq","freq_band","con_method"],
                                                output_names = "conmat_file",
                                                function = spectral_proc),
        
from neuropype_ephy.spectral import spectral_proc

class SpectralConnInputSpec(BaseInterfaceInputSpec):
    
    ts_file = traits.List(File(exists=True), desc='nodes * time series in .npy format', mandatory=True)
    
    sfreq = traits.Float(desc='sampling frequency', mandatory=True)
    
    freq_band = traits.List(traits.Float(exists=True), desc='frequency bands', mandatory=True)
    
    con_method = traits.Enum("coh","imcoh","plv","pli",desc='metric computed on time series for connectivity')
    
    
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
        
        self.conmat_file = spectral_proc(ts_file,sfreq,freq_band,con_method):

        return runtime
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["conmat_file"] = conmat_file
        
        return outputs
        
