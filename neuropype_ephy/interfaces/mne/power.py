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
    epochs_file = traits.File(exists=True, desc='File with mne.Epochs', mandatory=True)
    fmin = traits.Float(desc='lower psd frequency', mandatory=False)
    fmax = traits.Float(desc='higher psd frequency', mandatory=False)
    method = traits.Enum('welch', 'multitaper', desc='power spectral density computation method')


class PowerOutputSpec(TraitedSpec):
    psds_file = File(exists=True, desc='psd tensor in .npy format')
    freqs_file = File(exists=True, desc='freqs array in .npy format')

class Power(BaseInterface):
    """
    Compute power spectral density on epochs
    """
    input_spec = PowerInputSpec
    output_spec = PowerOutputSpec

    def _run_interface(self, runtime):
        print 'in Power'
        epochs_file = self.inputs.epochs_file
        fmin = self.inputs.fmin
        fmax = self.inputs.fmax
        method = self.inputs.method
        self.psds_file, self.freqs_file = compute_and_save_psd(epochs_file, fmin, fmax, method)
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['psds_file'] = self.psds_file
        outputs['freqs_file'] = self.freqs_file
        return outputs