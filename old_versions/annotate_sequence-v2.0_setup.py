from setuptools import setup, find_packages

import os 

setup(
    name="annotate_sequence",
    version="2.0",
    packages=["annotate_sequence-v2.0"],
    author="James Boocock and Tanya Flynn",
    author_email="james.boocock@otago.ac.nz",
    description="Annotate a sequence with variants from a .vcf file",
    license="Mit",
    zip_safe=False,
     entry_points={
        'console_scripts': [
            'annotate_sequence-v2.0 = annotate_sequence.loc_filter:main',
        ]
        },
    url="github.com/theboocock/snp_design",
    use_2to3=True,
)
Status API Training Shop Blog About
