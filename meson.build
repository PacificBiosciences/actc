project(
  'actc',
  'cpp',
  version : '0.0.2',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++17',
    'b_ndebug=if-release'],
  license : 'BSD-3',
  meson_version : '>= 0.52.0')

# CXXFLAGS
actc_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
  '-Werror=strict-aliasing',
  '-Wno-unused-parameter',
  '-Wno-unused-local-typedef']
    if cpp.has_argument(cflag)
      actc_flags += cflag
    endif
endforeach

# dependencies #
## threads
actc_thread_dep = dependency('threads', required : true)
## boost
actc_boost_dep = dependency('boost', required : true)
## zlib
actc_zlib_dep = dependency('zlib', required : true, fallback : ['zlib', 'zlib_dep'])
## pbbam
actc_pbbam_dep = dependency('pbbam', fallback : ['pbbam', 'pbbam_dep'])
## pbcopper
actc_pbcopper_dep = dependency('pbcopper', fallback : ['pbcopper', 'pbcopper_dep'])
# htslib
actc_htslib_dep = dependency('htslib', required : true, version : '>=1.4', fallback : ['htslib', 'htslib_dep'])
# pancake
actc_pancake_dep = dependency('pancake', fallback : ['pancake', 'pancake_dep'])

actc_lib_deps = [
  actc_thread_dep,
  actc_boost_dep,
  actc_zlib_dep,
  actc_pbbam_dep,
  actc_pbcopper_dep,
  actc_htslib_dep,
  actc_pancake_dep,
]

actc_src_include_directories = include_directories('src')

# sources + executable
actc_main = executable(
  'actc',
  files([
    'src/main.cpp',
    'src/AlignmentResult.cpp',
    'src/AlignerUtils.cpp',
    'src/PancakeAligner.cpp',
  ]),
  install : true,
  dependencies : actc_lib_deps,
  include_directories : actc_src_include_directories,
  cpp_args : actc_flags)
