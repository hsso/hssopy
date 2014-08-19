#!/usr/bin/env python

from setuptools import setup

setup(name='hsso',
    description='Herschel Data Analysis Utilities',
    author='Miguel de Val-Borro',
    author_email='valborro@princeton.edu',
    version='0.1',
    packages=['hsso'],
    requires=[
        "numpy",
        "matplotlib",
    ],
    )
