"""
Docstring here
"""


class PreprocIputSpec(BaseInterfaceInputSpec):

class PreprocOutputSpec(TraitedSpec):


class Preproc(BaseInterface):
    """
    Perform preprocessing and ICA artifacts removal
    """
    def _run_interface(self, runtime):
        return runtime

    def _list_outputs(self):
        return outputs