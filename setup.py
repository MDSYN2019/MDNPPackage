#!/usr/bin/env python3

"""
Author: Sang Young Noh
----------------------
Last Updated: 02/07/2022
------------------------
"""

"""
The setup.py file is at the heart of a python project. It describes 
all of the metadata about your project. There are quite a few fields
you can add to a project to give it a rich set of metadata describing
the project.

However, there are only three required names: 

- name 
- version
- packages 

The name field must be unique if you wish to publish your pakcage 

"""
import os
from setuptools import setup
from setuptools import find_packages

# get key package details from py_pkg/__version__.py
about = {}  # type: ignore
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'py_pkg', '__version__.py')) as f:
    exec(f.read(), about)

# load the README file and use it as the long_description for PyPI
with open('README.md', 'r') as f:
    readme = f.read()

# package configuration - for reference see:
# https://setuptools.readthedocs.io/en/latest/setuptools.html#id9

setup(
    name='Martini-PyNP', # name of the package 
    description='Python package to create simulation templates for \ 
    Nanoparticle mixed simulations', # Short description 
    long_description=open("README.txt").read(), # Long description
    long_description_content_type='text/markdown',
    version='0.1.1',
    author='Sang Young Noh',
    author_email='sangyoung123@googlemail.com',
    url=about['__url__'],
    packages=find_packages('MDNPPackage', exclude = ['test*']), # Do not include any folder that include tests*
    include_package_data=True,
    python_requires=">=3.7.*",
    install_requires=['numpy', 'requests', 'mdanalysis', 'mdtraj', 'vermouth', 'scipy', 'rdkit-pypi', 'ParmEd'],
    license='MIT',
    zip_safe=False, # ?

    setup_requires = ['pytest-runner'],

    tests_require = ['pytest==4.4.1'],

    test_suite = 'tests',

    entry_points={
        'console_scripts': ['py-package-template=py_pkg.entry_points:main'],
    },

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Programming Language :: Python :: 3.7',
    ],
    
    keywords='package development template'
)
