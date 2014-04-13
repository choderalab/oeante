"""oeante: OpenEye toolki-based Antechamber implementation for parameterizing small molecules with GAFF.

Based on GPLv2 version by Richard Dixon.
"""
from __future__ import print_function
import os
import sys
from os.path import relpath, join
from setuptools import setup, find_packages
DOCLINES = __doc__.split("\n")

########################
__version__ = '1.0'
VERSION = __version__
ISRELEASED = False
########################
CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GPLv2
Programming Language :: Python
Topic :: Scientific/Engineering :: Bio-Informatics
Topic :: Scientific/Engineering :: Chemistry
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""


def find_package_data():
    files = []
    for root, dirnames, filenames in os.walk('oeante'):
        for fn in filenames:
            files.append(relpath(join(root, fn), 'oeante'))

    return files

setup(
    name='oeante',
    author='Richard Dixon and John Chodera',
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    version=__version__,
    license='GPLv2',
    url='https://github.com/choderalab/oeante',
    platforms=['Linux', 'Mac OS-X', 'Unix'],
    classifiers=CLASSIFIERS.splitlines(),
    packages=find_packages(),
    package_data={'oeante': find_package_data()},
    zip_safe=False,
    install_requires=[],
    entry_points={'console_scripts': ['oeante = oeante.oeante:main']})

