"""
Power Pipeline
"""
# Authors: Dmitrii Altukhov <dm-altukhov@ya.ru>
#          Annalisa Pascarella <a.pascarella@iac.cnr.it>

import nipype.pipeline.engine as pe

from nipype.interfaces.utility import IdentityInterface
from neuropype_ephy.interfaces.mne.power import Power


def create_pipeline_power(main_path, pipeline_name='power',
                          fmin=0, fmax=300, method='welch',
                          is_epoched=False):
    """
    Description:
    
        Wraps functions of MNE to compute PSD of epoch or raw data
    """    

    pipeline = pe.Workflow(name=pipeline_name)
    pipeline.base_dir = main_path

    print '*** main_path -> %s' % main_path + ' ***'

    # define the inputs of the pipeline
    inputnode = pe.Node(IdentityInterface(fields=['fif_file']),
                        name='inputnode')

    power_node = pe.Node(interface=Power(), name='power')
    power_node.inputs.fmin = fmin
    power_node.inputs.fmax = fmax
    power_node.inputs.method = method
    power_node.inputs.is_epoched = is_epoched

    pipeline.connect(inputnode, 'fif_file', power_node, 'epochs_file')

    return pipeline
