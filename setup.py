#!/usr/bin/env python

from setuptools import setup
from peako import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name="peako",
      version=__version__,
      description="A package for discovering motifs in ChIP-seq datasets with knockout controls",
      author="Danielle Denisko",
      author_email="danielle.denisko@mail.utoronto.ca",
      license='GPLv3',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url="https://bitbucket.org/hoffmanlab/peako",
      packages=['peako'],
      package_data={'peako': ['data/*']},
      entry_points = {
        'console_scripts': [
            'peako=peako.peako_snakemake_wrapper:main',
            ],
        },
      zip_safe=False)
