
from setuptools import setup, find_packages
from setuptools import Extension
import glob

setup(
    name='spirecon',
    version='0.1.0',
    author='Yigit Dallilar',
    author_email='ydallilar@mpe.mpg.de',
    packages = find_packages(),
    url='',
    description='Unofficial SPIFFIER Reconstruction',
    package_data={'spirecon' : ['calib/*.txt']},
    scripts=glob.glob('scripts/*.py'),
    python_requires="~=3.6",
    install_requires=['numpy', 'scipy', 'astropy', 'matplotlib'],
)