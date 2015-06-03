from setuptools import setup, find_packages

#here = path.abspath(path.dirname(__file__))
# Get the long description from the relevant file
#with open(path.join(here, 'DESCRIPTION.rst'), encoding='utf-8') as f:
#long_description = f.read()

setup(
        name = "NaviCom",
        version = "0.6",
        description = "Package for data visualisation in NaviCell",
        long_description = "Package for data visualisation in NaviCell",
        url = "", #TODO
        author = "Mathurin Dorel",
        author_email = "mathurin.dorel@curie.fr",
        license = "GPL", # TODO See with Andrei
        classifiers = [
            "Development Status :: 3 - Alpha",
            "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Programming Language :: Python :: 3"
            ],
        keywords = "navicell high-throughput_data",
        packages = find_packages(),
        install_requires = [
            "numpy>=1.4.1",
            "re>=2.2.1"
            ],
        data_files = []
    )
