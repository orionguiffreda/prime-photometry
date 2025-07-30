from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name='photometrus',
    version='0.0.1',
    description='Photometrus is the photometry pipeline for the PRIME telescope',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/orionguiffreda/prime-photometry',
    author='Orion Guiffreda',
    author_email='oriogui@umd.edu',
    keywords='photometry',
    packages=find_packages('photometrus'),
)
