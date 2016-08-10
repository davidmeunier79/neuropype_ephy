"""
Power computation module
"""

from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec
from nipype.utils.filemanip import split_filename

import nibabel as nbconvert
import numpy as np
import os

from neuropype_ephy.power import compute_and_save_psd

class PowerInputSpec(BaseInterfaceInputSpec):
    pass
    # epochs_file = traits.File(exists=True, desc='', mandatory=True)
    # fmin = traits.Float(desc='lower psd frequency')
    # mfmax = traits.Float(desc='higher psd frequency')



class PowerOutputSpec(TraitedSpec):
    pass
    # npy_file = File(exists=True, desc='psd tensor in .npy format')

class Power(BaseInterface):
    """
    Perform preprocessing and ICA artifacts removal
    """
    input_spec = PowerInputSpec
    output_spec = PowerOutputSpec

    def _run_interface(self, runtime):

        print 'in Power'
        compute_and_save_psd()
        return runtime

    def _list_outputs(self):
        return