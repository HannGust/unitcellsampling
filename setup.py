#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages, Extension
from setuptools_rust import Binding, RustExtension

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    "ase==3.20.1",
    "mendeleev==0.2.17",
    "pandas==1.1.2",
    "sqlalchemy==1.3.19",
    "scikit-learn==0.23.2",
    "python-dotenv==0.14.0",
    "pathlib==1.0.1",
    "numpy==1.19.2",
    "setuptools-rust==0.11.5",
    "pymatgen==2020.12.3",
    "psutil==5.8.0",
    "gemmi==0.4.4",
    "persist-queue==0.6.0"
]

setup_requirements = []

test_requirements = []

setup(
    author="Benjamin Bolbrinker",
    author_email='benjamin.bolbrinker@kemi.uu.se',
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Package for efficient calculation and perturbation of 3D energy grids.",
    entry_points={
        'console_scripts': [
            'unitcellsampling=unitcellsampling.cli:main',
            'grid-correct-mmc=unitcellsampling.cli:run_mmc',
            'grid-average-by-symmetry=unitcellsampling.cli:average_grid_values_by_symmetry',
            'grid-downsample=unitcellsampling.cli:down_sample',
            'grid-get-partial-charge=unitcellsampling.cli:get_partial_ion_charges',
            'grid-correct-mmc-create-traj=unitcellsampling.cli:view_trajectory',
        ],
    },
    install_requires=requirements,
    license="GNU General Public License v3",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='unitcellsampling',
    name='unitcellsampling',
    packages=find_packages(include=['unitcellsampling', 'unitcellsampling.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/benjaminbolbrinker/unitcellsampling',
    version='0.1.0',
    zip_safe=False,
    rust_extensions=[RustExtension(
        "gridmc", path="src/gridmc/Cargo.toml", binding=Binding.PyO3)]
)
