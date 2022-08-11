### Build the cython extension ###

# compile .c for fastc.pyx
fastc = custom_target('fastc',
  output : 'fastc.c',
  input : 'fastc.pyx',
  command : [cython, '-3', '--fast-fail', '@INPUT@', '-o', '@OUTPUT@']
)

npymath_path = incdir_numpy / '..' / 'lib'
npymath_lib = meson.get_compiler('c').find_library('npymath', dirs: npymath_path)
npyrandom_path = incdir_numpy / '..' / '..' / 'random' / 'lib'
npyrandom_lib = meson.get_compiler('c').find_library('npyrandom', dirs: npyrandom_path)

# generate .so
py3.extension_module(
  'fastc', fastc,
  include_directories: [inc_np, npymath_path, npyrandom_path],
  dependencies : [py3_dep, npymath_lib, npyrandom_lib],
  install: true,
  install_dir: 'LaNeta'
)
