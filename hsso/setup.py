#!/usr/bin/env python

from distutils.core import setup

setup(name='hsso',
    description='Herschel Data Analysis Utilities',
    author='Miguel de Val-Borro',
    author_email='valborro@princeton.edu',
    version='0.1',
    py_modules=['class_utils'],
    install_requires=[
        "numpy",
        "matplotlib",
    ],
    )
