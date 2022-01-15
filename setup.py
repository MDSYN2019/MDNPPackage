#!/usr/bin/env python3

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
    name='Martini-PyNP', # 
    description='Python package to create simulation templates for Nanoparticle mixed simulations',
    long_description=readme,
    long_description_content_type='text/markdown',
    version='0.0.1',
    author='Sang Young Noh',
    author_email='sangyoung123@googlemail.com',
    url=about['__url__'],
    packages=['mdanalysis'],
    include_package_data=True,
    python_requires=">=3.7.*",
    install_requires=['numpy', 'requests'],
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
