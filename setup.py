from setuptools import setup, find_packages
import pylastic

VERSION = pylastic.__version__
DESCRIPTION = 'Tensor Analysis Library'
LONG_DESCRIPTION = 'Library for analyzing the crystal stiffness tensor'
with open('requirements.txt', 'r') as fio:
    REQUIREMENTS = fio.read().split('\n')
# Setting up
setup(
    # the name must match the folder name 'verysimplemodule'
    name="pylastic",
    version=VERSION,
    author="Roman Fakhrutdinov",
    author_email="<summedjesters@gmail.com>",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=REQUIREMENTS,

)