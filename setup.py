from setuptools import setup, find_packages

setup(
    name='pyrrhenius',
    version='1.0.0',
    description='A python package for calculating mineral electric conductivity',
    author='Kevin A. Mendoza',
    author_email='kevinmendoza@icloud.com',
    packages=find_packages(exclude=['tests','publication_images',
                                    'mineral_ensembles','ensemble_images','doc']),
    install_requires=[
        'numpy',
        'pandas',
        'pyfluids',
        'scipy'
    ],
)