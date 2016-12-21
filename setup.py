#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup


setup(
    name='neuropype_ephy',
    version='0.0.1dev',
    packages=['neuropype_ephy'],
    author=['David Meunier',
            'Annalisa Pascarella',
            'Dmitrii Altukhov'],
    description='Definition of functions used\
                 as Node for electrophy (EEG/MEG)\
                 pipelines within nipype framework',
    lisence='MIT',
    install_requires=['numpy>=1.3.0', 'mne']
)
