# -*- coding: utf-8 -*-
"""
All nodes for import that are NOT specific to a ephy package
"""
import numpy as np
import os

from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec, isdefined
    
from nipype.utils.filemanip import split_filename as split_f
    
############################################################################################### ImportMat #####################################################################################################

from neuropype_ephy.import_mat import import_tsmat_to_ts

class ImportMatInputSpec(BaseInterfaceInputSpec):
    
    tsmat_file = traits.File(exists=True, desc='nodes * time series in .mat (matlab format format', mandatory=True)
    
    data_field_name = traits.String("F", desc='Name of the structure in matlab', usedefault=True)
    
    good_channels_field_name = traits.String('ChannelFlag', desc='Boolean structure for choosing nodes, name of structure in matlab file', mandatory=False)
    
class ImportMatOutputSpec(TraitedSpec):
    
    ts_file = File(exists=True, desc="time series in .npy format")
    
class ImportMat(BaseInterface):
    
    """
    Compute spectral connectivity in a given frequency bands
    """
    input_spec = ImportMatInputSpec
    output_spec = ImportMatOutputSpec

    def _run_interface(self, runtime):
                
        print 'in ImportMat'
        
        tsmat_file = self.inputs.tsmat_file
        
        data_field_name = self.inputs.data_field_name
        
        good_channels_field_name = self.inputs.good_channels_field_name

        if not isdefined(good_channels_field_name):
            good_channels_field_name = None
            
        self.ts_file = import_tsmat_to_ts(ts_file,data_field_name,good_channels_field_name)
        
        return runtime
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["ts_file"] = self.ts_file
        
        return outputs
        