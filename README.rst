.. _readme:

README
******

Description
===========

Neuropype package of functions for electrophysiology analysis, can be used from
neuropype_graph and nipype


Documentation
=============

https://annapasca.github.io/neuropype/ephypype/neuropype_ephy.html


Installation
=============

Requirements
------------

Up to now neuropype_ephy works only with python2; python3 compatibility is planned for later releases

* numpy
* scikit-learn
* mne
* nipype
* neuropype_graph

Some of these dependencies (numpy, scikit-learn) you should install manually, others (mne, nipype) 
are installed automatically during neuropype_ephy installation (see "Install neuropype_ephy").
To install neuropype_graph see "Install neuropype_graph"

Install package
---------------

Install neuropype_ephy
++++++++++++++++++++++

::

    git clone https://github.com/annapasca/neuropype_ephy.git
    cd neuropype_ephy
    sudo python setup.py develop
    cd ..


Install neuropype_graph
+++++++++++++++++++++++

:: 

    git clone https://github.com/davidmeunier79/neuropype_graph.git
    cd neuropype_graph
    pip install .
    cd ..

see |README_graph| for more information.

.. |README_graph| raw:: html

   <a href="http://davidmeunier79.github.io/neuropype_graph/includeme.html" target="_blank">README</a>

Software
--------

Freesurfer
++++++++++
1. Download Freesurfer sotware:

https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall

2. Follow the Installation instructions

https://surfer.nmr.mgh.harvard.edu/fswiki/LinuxInstall


MNE
+++

1. Download MNE sotware:

http://martinos.org/mne/dev/install_mne_c.html

