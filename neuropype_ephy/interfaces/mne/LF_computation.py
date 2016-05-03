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


class LFComputationConnInputSpec(BaseInterfaceInputSpec):

    sbj_id = traits.String(desc='subject id', mandatory=True)

    sbj_dir = traits.Directory(exists=True, desc='Freesurfer main directory',
                               mandatory=True)

    raw_info = traits.Any(desc='raw info', mandatory=True)

    spacing = traits.String('ico-5', desc='spacing to use to setup a source space',
                            mandatory=False)

    aseg = traits.Bool(False, desc='if true sub structures will be considered',
                       mandatory=False)

    aseg_labels = traits.List(desc='list of substructures in the src space',
                              mandatory=False)


class LFComputationConnOutputSpec(TraitedSpec):

    fwd_filename = File(exists=False, desc='LF matrix')


class LFComputation(BaseInterface):
    """
    Compute LF matrix using MNE Python functions
    """
    input_spec = LFComputationConnInputSpec
    output_spec = LFComputationConnOutputSpec

    def _get_fwd_filename(self, raw_info, aseg, spacing):

        data_path, raw_fname, ext = split_f(raw_info['filename'])

        if aseg == traits.Undefined:
            fwd_filename = op.join(data_path, '%s-%s-fwd.fif'
                                   % (raw_fname, spacing))
        else:
            fwd_filename = op.join(data_path, '%s-%s-aseg-fwd.fif'
                                   % (raw_fname, spacing))

        print '*** fwd_filename %s ***' % fwd_filename
        return fwd_filename

    def _run_interface(self, runtime):

        sbj_id = self.inputs.sbj_id
        sbj_dir = self.inputs.sbj_dir
        raw_info = self.inputs.raw_info
        aseg = self.inputs.aseg
        spacing = self.inputs.spacing
        aseg_labels = self.inputs.aseg_labels

        if aseg == traits.Undefined:
            self.fwd_filename = self._get_fwd_filename(raw_info, aseg, spacing)
        else:
            self.fwd_filename = self._get_fwd_filename(raw_info, aseg, spacing)

        # check if we have just created the fwd matrix
        if not op.isfile(self.fwd_filename):
            bem = create_bem_sol(sbj_dir, sbj_id)  # bem solution

            src = create_src_space(sbj_dir, sbj_id, spacing)  # src space

            if aseg != traits.Undefined:
                src = create_mixed_source_space(sbj_dir, sbj_id, spacing,
                                                aseg_labels, src)

            n = sum(src[i]['nuse'] for i in range(len(src)))
            print('il src space contiene %d spaces e %d vertici'
                  % (len(src), n))

            trans_fname = is_trans(raw_info)

            # TODO: ha senso una funzione con un solo cmd?
            compute_fwd_sol(raw_info, trans_fname, src, bem, self.fwd_filename)
        else:
            print '*** FWD file %s exists!!!' % self.fwd_filename

        return runtime

    def _list_outputs(self):

        outputs = self._outputs().get()

        outputs['fwd_filename'] = self.fwd_filename

        return outputs
