import setuptools

#with open('README.md', 'r') as fh:
#   long_description = fh.read()

setuptools.setup(
    name='SRAscraper',
    version='0.1.1',
    author='Jake Lehle',
    author_email='jlehle@txbiomed.org',
    description='Snakmake pipeline to pull fastq sequencing files from the SRA database using their GSE accession identifiers from GEO.',
    packages=[
        'cli',
        'snakemake_wrapper'
    ],
    py_modules=[
        'cli'
    ],
    package_data={
        'snakemake_wrapper': [
            'scripts/*.py',
            'scripts/*.sh',
            'envs/*.yaml',
            'Snakefile'
        ],
        'cli': [
            'optionals.yaml'
        ]
    },
    install_requires=[
        'click',
        'ruamel.yaml',
        'snakemake'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research'
    ],
    entry_points='''
        [console_scripts]
        SRAscraper=cli.cli:main
    '''
)
