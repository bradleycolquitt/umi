from distutils.core import setup
from Cython.Distutils import build_ext
from Cython.Build import cythonize

setup(
    name = 'Barcode splitter',
    ext_modules = cythonize("*.pyx", gdb_debug=True)
)
