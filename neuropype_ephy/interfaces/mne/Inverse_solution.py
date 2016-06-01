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

from neuropype_ephy.compute_inv_problem import compute_ROIs_inv_sol
from neuropype_ephy.preproc import create_reject_dict

from mne import compute_raw_covariance, pick_types, write_cov
from mne.io import Raw

class InverseSolutionConnInputSpec(BaseInterfaceInputSpec):

    sbj_id = traits.String(desc='subject id', mandatory=True)

    sbj_dir = traits.Directory(exists=True, desc='Freesurfer main directory',
                               mandatory=True)

    raw_filename = traits.File(exists=True, desc='raw filename', mandatory=True)

    cov_filename = traits.File(exists=True, desc='Noise Covariance matrix',
                               mandatory=True)

    fwd_filename = traits.File(exists=True, desc='LF matrix', mandatory=True)

    is_epoched = traits.Bool(desc='if true raw data will be epoched',
                             mandatory=False)    
                             
    event_id = traits.Int(None, desc='the id of the event to consider.', mandatory=False)
    
    t_min = traits.Float(None, desc='start time before event', mandatory=False)

    t_max = traits.Float(None, desc='end time after event', mandatory=False)
                   
    inv_method = traits.String(desc='possible inverse methods are \
                               sLORETA, MNE, dSPM', mandatory=True)

    snr = traits.Float(1.0, usedefault=True, desc='use smaller SNR for \
                       raw data', mandatory=False)

    parc = traits.String('aparc', usedefault=True,
                         desc='the parcellation to use: aparc vs aparc.a2009s',
                         mandatory=False)

    aseg = traits.Bool(desc='if true sub structures will be considered',
                       mandatory=False)

    aseg_labels = traits.List(desc='list of substructures in the src space',
                              mandatory=False)

    is_blind = traits.Bool(desc='if in the source space there are ROI removed',
                           mandatory=False)

    labels_removed = traits.List(desc='list of label we consider in the blind case',
                                 mandatory=False)


class InverseSolutionConnOutputSpec(TraitedSpec):

    ts_file = File(exists=False, desc='source reconstruction in .npy format')
    labels = File(exists=False, desc='labels file in pickle format')
    label_names = File(exists=False, desc='labels name file in txt format')
    label_coords = File(exists=False, desc='labels coords file in txt format')


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
        raw_filename = self.inputs.raw_filename
        cov_filename = self.inputs.cov_filename
        fwd_filename = self.inputs.fwd_filename
        is_epoched = self.inputs.is_epoched
        event_id = self.inputs.event_id
        t_min = self.inputs.t_min
        t_max = self.inputs.t_max
        inv_method = self.inputs.inv_method
        snr = self.inputs.snr
        parc = self.inputs.parc
        aseg = self.inputs.aseg
        aseg_labels = self.inputs.aseg_labels
        is_blind = self.inputs.is_blind
        labels_removed = self.inputs.labels_removed

        self.ts_file, self.labels , self.label_names, self.label_coords= compute_ROIs_inv_sol(raw_filename, sbj_id, sbj_dir,
                                                                                              fwd_filename,
                                                                                              cov_filename,
                                                                                              is_epoched,
                                                                                              event_id, t_min, t_max,
                                                                                              snr, inv_method, parc,
                                                                                              aseg, aseg_labels,
                                                                                              is_blind, labels_removed)

        return runtime

    def _list_outputs(self):

        outputs = self._outputs().get()

        outputs['ts_file'] = self.ts_file
        outputs['labels'] = self.labels
        outputs['label_names'] = self.label_names
        outputs['label_coords'] = self.label_coords

        return outputs


class NoiseCovarianceConnInputSpec(BaseInterfaceInputSpec):

    cov_fname_in = traits.File(exists=False, desc='file name for Noise Covariance Matrix')

    raw_filename = traits.File(exists=True, desc='raw data filename')


class NoiseCovarianceConnOutputSpec(TraitedSpec):

    cov_fname_out = File(exists=False, desc='LF matrix')


class NoiseCovariance(BaseInterface):
    """
    Compute the inverse solution on raw data considering N_r regions in source
    space based on a FreeSurfer cortical parcellation
    """
    input_spec = NoiseCovarianceConnInputSpec
    output_spec = NoiseCovarianceConnOutputSpec

    def _run_interface(self, runtime):

        raw_filename = self.inputs.raw_filename
        cov_fname_in = self.inputs.cov_fname_in
        
        if cov_fname_in == '' or not op.exists(cov_fname_in):
            
            raw = Raw(raw_filename)
            
            data_path, basename, ext = split_f(raw.info['filename'])
            self.cov_fname_out = op.join(data_path, '%s-cov.fif' % basename)

            print '\n*** COMPUTE RAW COV ***\n' + self.cov_fname_out

#            reject = dict(mag=4e-12, grad=4000e-13, eog=250e-6)
            reject = create_reject_dict(raw.info)

            picks = pick_types(raw.info, meg=True, ref_meg=False,
                               exclude='bads')

            # TODO: which method? method='auto' too long for raw data
            noise_cov = compute_raw_covariance(raw, picks=picks, reject=reject)

            write_cov(self.cov_fname_out, noise_cov)

        else:
            print '\n *** NOISE cov file %s exists!!! \n' % cov_fname_in
            self.cov_fname_out = cov_fname_in

        return runtime

    def _list_outputs(self):

        outputs = self._outputs().get()

        outputs['cov_fname_out'] = self.cov_fname_out

        return outputs
