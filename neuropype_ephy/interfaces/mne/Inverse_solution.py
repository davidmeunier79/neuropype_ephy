# -*- coding: utf-8 -*-
"""
Created on Mon May  2 21:30:03 2016

@author: karim
"""

# -*- coding: utf-8 -*-
"""
Created on Mon May  2 17:24:00 2016

@author: pasca
"""

# -*- coding: utf-8 -*-
import os.path as op

from nipype.utils.filemanip import split_filename as split_f

from nipype.interfaces.base import BaseInterface, BaseInterfaceInputSpec
from nipype.interfaces.base import traits, File, TraitedSpec

from neuropype_ephy.compute_fwd_problem import create_mixed_source_space
from neuropype_ephy.compute_fwd_problem import create_bem_sol, create_src_space
from neuropype_ephy.compute_fwd_problem import is_trans, compute_fwd_sol


class InverseSolutionConnInputSpec(BaseInterfaceInputSpec):

    sbj_id = traits.String(desc='subject id', mandatory=True)

    sbj_dir = traits.Directory(exists=True, desc='Freesurfer main directory',
                               mandatory=True)

    raw = traits.Any(desc='raw data', mandatory=True)

    snr = traits.Float(desc='use smaller SNR for raw data',
                       mandatory=False)

    parc = traits.String(desc='the parcellation to use: aparc vs aparc.a2009s',
                         mandatory=True)

    aseg_labels = traits.List(desc='list of substructures in the src space',
                              mandatory=False)

    fwd_filename = traits.File(exists=True, desc='LF matrix', mandatory=True)

    inv_method = traits.String(desc='inverse method used', mandatory=True)


class InverseSolutionConnOutputSpec(TraitedSpec):

    ts_file = File(exists=False, desc='source reconstruction in .npy format')
    labels = traits.Any(desc='labels list')


class InverseSolution(BaseInterface):
    """
    Compute the inverse solution on raw data considering N_r regions in source
    space based on a FreeSurfer cortical parcellation
    """
    input_spec = InverseSolutionConnInputSpec
    output_spec = InverseSolutionConnOutputSpec

    def _run_interface(self, runtime):

        sbj_id = self.inputs.sbj_id
        sbj_dir = self.inputs.sbj_dir
        raw = self.inputs.raw
        snr = self.inputs.snr
        parc = self.inputs.parc
        aseg_labels = self.inputs.aseg_labels
        fwd_filename = self.inputs.fwd_filename
        inv_method = self.inputs.inv_method



        return runtime

    def _list_outputs(self):

        outputs = self._outputs().get()

        outputs['ts_file'] = self.ts_file
        outputs['labels'] = self.labels
        
        return outputs