"""Interfaces for preprocessing nodes"""

from nipype.interfaces.base import BaseInterface,\
    BaseInterfaceInputSpec, traits, TraitedSpec

from neuropype_ephy.preproc import compute_ica


class CompIcaInputSpec(BaseInterfaceInputSpec):
    """Input specification for CompIca"""
    fif_file = traits.File(exists=True,
                           desc='raw meg data in fif format',
                           mandatory=True)
    ecg_ch_name = traits.String(desc='name of ecg channel')
    eog_ch_name = traits.String(desc='name of eog channel')
    n_components = traits.Float(desc='number of ica components')


class CompIcaOutputSpec(TraitedSpec):
    """Output specification for CompIca"""
    ica_file = traits.File(exists=True,
                           desc='file with ica solution in .fif',
                           mandatory=True)

    ica_ts_file = traits.File(exists=True,
                              desc='file with ica components in .fif',
                              mandatory=True)
    report_file = traits.File(exists=True,
                              desc='ica report in .html',
                              mandatory=True)


class CompIca(BaseInterface):
    """Compute ICA solution on raw fif data"""
    input_spec = CompIcaInputSpec
    output_spec = CompIcaOutputSpec

    def _run_interface(self, runtime):
        fif_file = self.inputs.fif_file
        ecg_ch_name = self.inputs.ecg_ch_name
        eog_ch_name = self.inputs.eog_ch_name
        n_components = self.inputs.n_components

        ica_file, ica_ts_file, report_file = compute_ica(fif_file,
                                                         ecg_ch_name,
                                                         eog_ch_name,
                                                         n_components)
        self.ica_file = ica_file
        self.ica_ts_file = ica_ts_file
        self.report_file = report_file

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['ica_file'] = self.ica_file
        outputs['ica_ts_file'] = self.ica_ts_file
        outputs['report_file'] = self.report_file
        return outputs
