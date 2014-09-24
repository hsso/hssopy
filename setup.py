#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='hsso',
    description='Herschel Data Analysis Utilities',
    author='Miguel de Val-Borro',
    author_email='valborro@princeton.edu',
    version='0.1',
    packages=find_packages(),
    py_modules=['herschel'],
    scripts=['hsso/hipe.py',
        'wcslint.py'],
    requires=[
        "numpy",
        "matplotlib",
    ],
    )
