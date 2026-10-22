#from numpy.distutils.core import setup, Extension

from distutils.core import setup
from distutils.extension import Extension
import numpy

module = Extension(
	'symdali',
	include_dirs=[numpy.get_include()],
       sources = ['symdali.c']
         
)

with open("README.md", encoding="utf-8") as f:
    long_description = f.read()


setup (name = 'symdali',
       version = '1.0.0',
       description=(
       "Python C extension for the SymDALI framework "
       "to approximate gravitational-wave likelihoods"
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    ext_modules=[module],
    install_requires=["numpy>=1.26,<2"],
    python_requires=">=3.11",
    license_files=["LICENSE"],)