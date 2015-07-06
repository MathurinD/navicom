from setuptools import setup, find_packages
import sys

#here = path.abspath(path.dirname(__file__))
# Get the long description from the relevant file
#with open(path.join(here, 'DESCRIPTION.rst'), encoding='utf-8') as f:
#long_description = f.read()

version = sys.version_info
if (not version[0] == 3):
    raise ValueError("navicom is only valid for python 3")

setup(
        name = "NaviCom",
        version = "1.0",
        description = "Package for data visualisation in NaviCell",
        long_description = "Package for data visualisation in NaviCell",
        url = "", #TODO
        author = "Mathurin Dorel",
        author_email = "mathurin.dorel@curie.fr",
        license = "LGPL",
        classifiers = [
            "Development Status :: 3 - Alpha",
            "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Programming Language :: Python :: 3"
            ],
        keywords = "navicell high-throughput_data",
        packages = find_packages(),
        install_requires = [
            "numpy>=1.4.1"
            ],
        data_files = []
    )
