## python setup.py build_ext --inplace

# If you get this error: clang: error: unsupported option '-fopenmp'
# run the following:
# brew install llvm libomp
# export CC=/usr/local/opt/llvm/bin/clang
# export CXX=/usr/local/opt/llvm/bin/clang++

#from distutils.core import setup
#from Cython.Build import cythonize
#import numpy as np

#setup(
#    ext_modules = cythonize("PRIME_reduction.pyx"),
#    include_dirs = [np.get_include()],
#    extra_compile_args=["-fopenmp"],
#    extra_link_args=["-fopenmp"]
#    )


from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension(
        "ramp_reduce",
        ["ramp_reduce.pyx"],
        include_dirs=[np.get_include()],
        extra_compile_args=["-std=c99", "-fopenmp"],
        extra_link_args=["-fopenmp"],
    )
]

setup(
    ext_modules=cythonize(extensions, language_level=3),
)


