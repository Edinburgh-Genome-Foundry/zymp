import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open('zymp/version.py').read()) # loads __version__

setup(
    name='zymp',
    version=__version__,
    author='Zulko',
    url='https://github.com/Edinburgh-Genome-Foundry/zymp',
    description='Design compact sequences with many enzyme restriction sites.',
    long_description=open('pypi-readme.rst').read(),
    license='MIT',
    keywords="DNA sequence design restriction site array",
    packages=find_packages(exclude='docs'),
    install_requires=['numpy', 'dnachisel', 'dna_features_viewer', 'biopython',
                      'dnacauldron', 'proglog'])
