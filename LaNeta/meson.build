### Build the cython extension ###

# compile .c for fastc.pyx
fastc = custom_target('fastc',
  output : 'fastc.c',
  input : 'fastc.pyx',
  command : [cython, '-3', '--fast-fail', '@INPUT@', '-o', '@OUTPUT@']
)

npymath_path = incdir_numpy / '..' / 'lib'
npymath_lib = meson.get_compiler('c').find_library('npymath', dirs: npymath_path)

# generate .so
py3.extension_module(
  'fastc', fastc,
  include_directories: [inc_np, npymath_path],
  dependencies : [py3_dep, npymath_lib],
  install: true,
  subdir: 'LaNeta'
)

### Install stuff: python files + the extension just built ###

#create a list of files which will be installed as source
python_sources = [
    'ThLd.py'
]

# install source files in current dir
py3.install_sources(
    python_sources,
    pure: false,
    subdir: 'LaNeta'
)
