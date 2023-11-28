from setuptools import setup, find_packages

DESCRIPTION = 'Tensor Analysis Library'
LONG_DESCRIPTION = 'Library for analyzing the crystal stiffness tensor'

# Setting up
setup(
    # the name must match the folder name 'verysimplemodule'
    name="pylastic",
    version='0.1.0',
    author="Roman Fakhrutdinov",
    author_email="<summedjesters@gmail.com>",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=['aiohttp~=3.8.5', 'pydantic~=2.4.2', 'Mako~=1.2.4', 'numpy~=1.26.0', 'scipy~=1.11.2',
                    'requests~=2.31.0', 'plotly~=5.18.0', 'setuptools'],
    include_package_data=True,
    package_data={'': ['*.html']}
)
