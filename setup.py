#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'Bio==1.7.1',
    'biopython==1.83',
    'cyvcf2==0.30.28',
    'datatable==1.1.0',
    'icecream==2.1.3',
    'kipoiseq==0.7.1',
    'loguru==0.7.2',
    'matplotlib==3.9.0',
    'numpy',
    'pandas==2.2.2',
    'pyfaidx==0.8.1.1',
    'pyliftover==0.4.1',
    'pysam==0.22.0',
    'scikit_learn==1.4.1.post1',
    'scipy==1.13.1',
    'seaborn==0.13.2'
]

test_requirements = [
    'coverage==4.5.4',
    'mypy',  # Not specified in requirements_dev.txt but needed for type checking
    'pytest>=3',
    'ruff==0.3.5',
    'pip==19.2.3',
    'bump2version==0.5.11',
    'wheel==0.33.6',
    'watchdog==0.9.0',
    'tox==3.14.0',
    'Sphinx==7.2.6',
    'twine==5.0.0'
]

setup(
    author="Rohit Ghosh",
    author_email='rohit.ghosh@yale.edu',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Analysis of single-cell Massively Parallel Reporter Assay (MPRA) data",
    entry_points={
        'console_scripts': [
            'worktools=worktools.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    package_data={
        'worktools': ['data/*.mtx'],
    },
    keywords='worktools',
    name='worktools',
    packages=find_packages(include=['worktools', 'worktools.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/rghosh670/worktools',
    version='0.1.0',
    zip_safe=False,
)
