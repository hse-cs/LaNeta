# Declare global project options:
# package name
# programming languages
project(
  'LaNeta',
  ['c', 'cython'],
  default_options : ['buildtype=debugoptimized'],
  version: '3.1',
  license: 'GPL-3.0 License',
  meson_version: '>= 0.54',
)

# Seek the backend
if meson.backend() != 'ninja'
  error('Ninja backend required')
endif

# Seek for Cython installation
cython = find_program('cython', required : true)
if not cython.found()
  error('Cython not found.')
endif


# Python default installton path and basic dependecies
py_mod = import('python')
py3 = py_mod.find_installation()
py3_dep = py3.dependency()

# Thx rgommers "https://github.com/rgommers/scipy/blob/meson/scipy/meson.build"
# NumPy include directory - needed in all submodules
incdir_numpy = run_command(py3,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()
# Declare numpy as dependency
inc_np = include_directories(incdir_numpy)

numpy_nodepr_api = '-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION'

subdir('LaNeta')
