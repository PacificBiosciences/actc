actc_src_include_directories = include_directories('.')

actc_config = configuration_data()
actc_config.set('ACTC_VERSION', meson.project_version())

actc_gen_headers = [
  vcs_tag(
    input : 'LibraryGitHash.hpp.in',
    output : 'LibraryGitHash.hpp',
    command : ['git', 'describe', '--always', '--dirty=*'],
    replace_string : '@ACTC_GIT_SHA1@'),
  configure_file(
    input : 'LibraryVersion.hpp.in',
    output : 'LibraryVersion.hpp',
    configuration : actc_config),
]

# sources + executable
actc_main = executable(
  'actc',
  files([
    'main.cpp',
    'AlignmentResult.cpp',
    'AlignerUtils.cpp',
    'LibraryInfo.cpp',
    'PancakeAligner.cpp',
    'io/BamZmwReader.cpp',
    'io/BamZmwReaderConfig.cpp',
  ]) + actc_gen_headers,
  install : true,
  dependencies : actc_lib_deps,
  include_directories : actc_src_include_directories,
  cpp_args : actc_flags)
