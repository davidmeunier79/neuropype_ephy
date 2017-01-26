# -*- coding: utf-8 -*-
"""
Created on Mon May  2 17:24:00 2016

@author: pasca
"""

# -*- coding: utf-8 -*-
import os.path as op
import sys
import glob

from nipype.utils.filemanip import split_filename as split_f

from nipype.interfaces.base import BaseInterface, BaseInterfaceInputSpec
from nipype.interfaces.base import traits, File, TraitedSpec

from neuropype_ephy.compute_inv_problem import compute_ROIs_inv_sol
from neuropype_ephy.preproc import create_reject_dict
from mne import find_events, compute_raw_covariance, compute_covariance
from mne import pick_types, write_cov, Epochs
from mne.io import read_raw_fif


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
    
    events_id = traits.Dict(None, desc='the id of all events to consider.', mandatory=False)                         
    
    event_id = traits.Int(None, desc='the id of the event to consider.',
                          usedefault=None, mandatory=False)
    
    t_min = traits.Float(None, desc='start time before event', mandatory=False)

    t_max = traits.Float(None, desc='end time after event', mandatory=False)
                   
    is_evoked = traits.Bool(desc='if true if we want to analyze evoked data',
                             mandatory=False)

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

    save_stc = traits.Bool(desc='if true save stc', mandatory=False)


class InverseSolutionConnOutputSpec(TraitedSpec):

    ts_file = File(exists=False, desc='source reconstruction in .npy format')
    labels = File(exists=False, desc='labels file in pickle format')
    label_names = File(exists=False, desc='labels name file in txt format')
    label_coords = File(exists=False, desc='labels coords file in txt format')


class InverseSolution(BaseInterface):
    """
    Compute the inverse solution on raw or epoch data considering N_r regions
    in source space based on a FreeSurfer cortical parcellation
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
        is_evoked = self.inputs.is_evoked
        events_id = self.inputs.events_id
        inv_method = self.inputs.inv_method
        snr = self.inputs.snr
        parc = self.inputs.parc
        aseg = self.inputs.aseg
        aseg_labels = self.inputs.aseg_labels
        is_blind = self.inputs.is_blind
        labels_removed = self.inputs.labels_removed
        save_stc = self.inputs.save_stc

        self.ts_file, self.labels, self.label_names, self.label_coords = \
            compute_ROIs_inv_sol(raw_filename, sbj_id, sbj_dir, fwd_filename,
                                 cov_filename, is_epoched, event_id,
                                 t_min, t_max, is_evoked, events_id,
                                 snr, inv_method, parc,
                                 aseg, aseg_labels, is_blind, labels_removed,
                                 save_stc)

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

    is_epoched = traits.Bool(desc='true if we want to epoch the data',
                             mandatory=False)

    is_evoked = traits.Bool(desc='true if we want to analyze evoked data',
                            mandatory=False)

    events_id = traits.Dict(None, desc='the id of all events to consider.', mandatory=False)

    t_min = traits.Float(None, desc='start time before event', mandatory=False)

    t_max = traits.Float(None, desc='end time after event', mandatory=False)


class NoiseCovarianceConnOutputSpec(TraitedSpec):

    cov_fname_out = File(exists=False, desc='Noise covariances matrix')


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
        is_epoched = self.inputs.is_epoched
        is_evoked = self.inputs.is_evoked
        events_id = self.inputs.events_id
        t_min = self.inputs.t_min
        t_max = self.inputs.t_max

        data_path, basename, ext = split_f(raw_filename)

        self.cov_fname_out = op.join(data_path, '%s-cov.fif' % basename)

        if not op.isfile(cov_fname_in):
            if is_epoched and is_evoked:
                raw = read_raw_fif(raw_filename)
                events = find_events(raw)

                if not op.isfile(self.cov_fname_out):
                    print '\n*** COMPUTE COV FROM EPOCHS ***\n' + self.cov_fname_out

                    reject = create_reject_dict(raw.info)
                    picks = pick_types(raw.info, meg=True, ref_meg=False,
                                       exclude='bads')

                    epochs = Epochs(raw, events, events_id, t_min, t_max,
                                    picks=picks, baseline=(None, 0),
                                    reject=reject)

                    # TODO method='auto'? too long!!!
                    noise_cov = compute_covariance(epochs, tmax=0,
                                                   method='diagonal_fixed')
                    write_cov(self.cov_fname_out, noise_cov)
                else:
                    print '\n *** NOISE cov file %s exists!!! \n' % self.cov_fname_out
            else:
                '\n *** RAW DATA \n'
                # TODO creare una matrice diagonale?
                for er_fname in glob.glob(op.join(data_path, cov_fname_in)):
                    print '\n found file name %s  \n' % er_fname

                try:
                    if er_fname.rfind('cov.fif') > -1:
                        print '\n *** NOISE cov file %s exists!! \n' % er_fname
                        self.cov_fname_out = er_fname
                    else:
                        er_raw = read_raw_fif(er_fname)
                        if not op.isfile(self.cov_fname_out):
                            reject = create_reject_dict(er_raw.info)
                            picks = pick_types(er_raw.info, meg=True,
                                               ref_meg=False, exclude='bads')
    
                            noise_cov = compute_raw_covariance(er_raw, picks=picks,
                                                               reject=reject)
                            write_cov(self.cov_fname_out, noise_cov)
                        else:
                            print '\n *** NOISE cov file %s exists!!! \n' % self.cov_fname_out
                except NameError:
                    sys.exit("No covariance matrix as input!")

        else:
            print '\n *** NOISE cov file %s exists!!! \n' % cov_fname_in
            self.cov_fname_out = cov_fname_in

        return runtime

    def _list_outputs(self):

        outputs = self._outputs().get()

        outputs['cov_fname_out'] = self.cov_fname_out

        return outputs
