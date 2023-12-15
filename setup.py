from setuptools import setup
from setuptools import find_packages

setup(
    name="edNEGmodel_UQSA",
    packages=find_packages(),
    package_data={'data': ['initial_values/*.dat']},
    version="1.0")
