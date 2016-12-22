#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


setup(
    name='neuropype_ephy',
    version='0.0.1dev',
    packages=find_packages(),
    author=['David Meunier',
            'Annalisa Pascarella',
            'Dmitrii Altukhov'],
    description='Definition of functions used\
                 as Node for electrophy (EEG/MEG)\
                 pipelines within nipype framework',
    lisence='MIT',
    install_requires=[ 'mne',
                      'nipype',
                      'configparser']
)
