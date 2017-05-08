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

Up to now neuropype_ephy works only with **python2**; python3 compatibility is planned for later releases

* numpy
* scikit-learn
* ipywidgets
* matplotlib
* pandas
* mne
* nipype
* neuropype_graph

Some of these dependencies you should install manually (see :ref:`conda_install`), others are installed automatically
during neuropype_ephy installation (see :ref:`ephy_install`).
To install neuropype_graph see :ref:`graph_install`. 

We also recommend to install the  development master version of |MNE install| (see :ref:`mne_install`).

.. |MNE install| raw:: html

   <a href="http://martinos.org/mne/dev/install_mne_python.html#check-your-installation" target="_blank">MNE python</a>

.. note:: If you have Anaconda it is possible to create an environment using python2 by the command
	``conda create -n py27 python=2.7 ipykernel``

.. warning:: We also recommend to use the nipype version 0.12
	``pip install nipype==0.12``
   
Install package
---------------

.. _ephy_install:

Install neuropype_ephy
++++++++++++++++++++++

.. code-block:: bash

    git clone https://github.com/annapasca/neuropype_ephy.git
    cd neuropype_ephy
    pip install .
    cd ..


.. _graph_install:

Install neuropype_graph
+++++++++++++++++++++++

.. code-block:: bash 

    git clone https://github.com/davidmeunier79/neuropype_graph.git
    cd neuropype_graph
    pip install .
    cd ..

see |README_graph| for more information.

.. |README_graph| raw:: html

   <a href="http://davidmeunier79.github.io/neuropype_graph/includeme.html" target="_blank">README</a>


.. _mne_install:
   
Install MNE python
++++++++++++++++++

.. code-block:: bash 

    git clone git://github.com/mne-tools/mne-python.git
    cd mne-python
    sudo python setup.py develop
    cd ..

see |MNE install| for more information.


.. _conda_install:
   
Install dependencies with conda
+++++++++++++++++++++++++++++++

.. code-block:: bash 

    conda install pandas
    conda install ipywidgets
    conda install matplotlib
    conda install 'pyqt<5'


Software
--------

Freesurfer
++++++++++
1. Download Freesurfer software:

https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall

2. Follow the Installation instructions

https://surfer.nmr.mgh.harvard.edu/fswiki/LinuxInstall


MNE
+++

1. Download MNE software:

http://martinos.org/mne/dev/install_mne_c.html

