# -*- coding: utf-8 -*-
"""

Description:

All nodes for import that are NOT specific to a ephy package
"""
import numpy as np
import os

from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, TraitedSpec, isdefined

from nipype.interfaces.base import File
   
from nipype.utils.filemanip import split_filename as split_f
    
############################################################################################### ImportMat #####################################################################################################

from neuropype_ephy.import_mat import import_tsmat_to_ts

class ImportMatInputSpec(BaseInterfaceInputSpec):
    
    tsmat_file = traits.File(exists=True, desc='nodes * time series in .mat (matlab format format', mandatory=True)
    
    data_field_name = traits.String("F", desc='Name of the structure in matlab', usedefault=True)
    
    good_channels_field_name = traits.String('ChannelFlag', desc='Boolean structure for choosing nodes, name of structure in matlab file')
    
    hdf5_mat = traits.Bool(False, desc='mat file is in hdf5 mormat (> 7.3)', usedefault  = True)
    
class ImportMatOutputSpec(TraitedSpec):
    
    ts_file = traits.File(exists=True, desc="time series in .npy format")
    
class ImportMat(BaseInterface):
    
    """
    Description:
    
    Import matlab file to numpy ndarry, and save it as numpy file .npy
    
    Inputs:
    
    tsmat_file: 
        type = File, exists=True, desc='nodes * time series in .mat (matlab format format', mandatory=True
        
    data_field_name 
        type = String, default = "F", desc='Name of the structure in matlab', usedefault=True
        
    good_channels_field_name
        type = String, default = 'ChannelFlag', desc='Boolean structure for choosing nodes, name of structure in matlab file'
    
    hdf5_mat
		type = Bool, default = False, desc='mat file is in hdf5 mormat (> 7.3)', usedefault  = True
		
    Outputs:

    ts_file 
        type = File, exists=True, desc="time series in .npy format"
        
    """
    input_spec = ImportMatInputSpec
    output_spec = ImportMatOutputSpec

    def _run_interface(self, runtime):
                
        print 'in ImportMat'
        
        tsmat_file = self.inputs.tsmat_file
        
        data_field_name = self.inputs.data_field_name
        
        hdf5_mat = self.inputs.hdf5_mat
        
        good_channels_field_name = self.inputs.good_channels_field_name

        if not isdefined(good_channels_field_name):
            good_channels_field_name = None
            
        self.ts_file = import_tsmat_to_ts(tsmat_file,data_field_name,good_channels_field_name, hdf5_mat = hdf5_mat)
        
        return runtime
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["ts_file"] = self.ts_file
        
        return outputs
       
############################################################################################### ImportBrainVisionAscii #####################################################################################################

from neuropype_ephy.import_txt import split_txt

class ImportBrainVisionAsciiInputSpec(BaseInterfaceInputSpec):
    
    txt_file = File(exists = True, desc='Ascii text file exported from BrainVision', mandatory=True)
    
    sample_size  = traits.Float(desc = "Size (number of time points) of all samples", mandatory = True)
                             
    sep_label_name = traits.String("", desc='Separator between electrode name (normally a capital letter) and contact numbers', usedefault=True)
    
    repair = traits.Bool(True, desc='Repair file if behaves strangely  (adding space sometimes...)', usedefault  = True)
    
    sep = traits.Str(";",desc = "Separator between time points",usedefault = True)
    
class ImportBrainVisionAsciiOutputSpec(TraitedSpec):

    print "in ImportBrainVisionAsciiOutputSpec"
    splitted_ts_file = traits.File(exists=True, desc="splitted time series in .npy format")
    
    elec_names_file  = traits.File(exists=True, desc="electrode names in txt format")
    
class ImportBrainVisionAscii(BaseInterface):
    
    """
    Description:
    
    Import IntraEEG Brain Vision (unsplitted) ascii time series txt file and return splitted time series in .npy format, as well as electrode names in txt format
    
    Inputs:
    
    txt_file 
        type = File, exists=True, desc='Ascii text file exported from BrainVision', mandatory=True
    
    sample_size  
        type = Int, desc = "Size (number of time points) of all samples", mandatory = True
                             
    sep_label_name
        type = String, default = "", desc='Separator between electrode name (normally a capital letter) and contact numbers', usedefault=True
        
    repair
        type = Bool, default = True, desc='Repair file if behaves strangely  (adding space sometimes...)', usedefault  = True
    
    sep
        type = String, default = ";","Separator between time points",usedefault = True)
    
    Outputs:

    splitted_ts_file
        type  = File, exists=True, desc="splitted time series in .npy format"
    
    elec_names_file
        type = File, exists=True, desc="electrode names in txt format"
    
        
    """
    input_spec = ImportBrainVisionAsciiInputSpec
    output_spec = ImportBrainVisionAsciiOutputSpec

    def _run_interface(self, runtime):
                
        print 'in ImportBrainVisionAscii'
        
        txt_file = self.inputs.txt_file
        
        sample_size  = self.inputs.sample_size
                                
        sep_label_name = self.inputs.sep_label_name
        
        repair = self.inputs.repair
        
        sep = self.inputs.sep
        
        split_txt(txt_file = txt_file, sample_size = sample_size, sep_label_name = sep_label_name, repair = repair, sep = sep)
        
        return runtime
        
    def _list_outputs(self):
        
        outputs = self._outputs().get()
        
        outputs["elec_names_file"] = os.path.abspath("correct_channel_names.txt")

        outputs["splitted_ts_file"] = os.path.abspath("splitted_ts.npy")
        
        return outputs
        
        
