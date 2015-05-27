#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
import os

install_requires = ['numpy>=1.3.0',]

setup(
    name = "neuropype_ephy",
    version = '0.0.1dev',
    packages = ['neuropype_ephy'],
    install_requires=install_requires,
    author = "David Meunier",
    description = "Definition of function used as Node for electrophy ( EEG/MEG) pipelines within nipype framework"
)

