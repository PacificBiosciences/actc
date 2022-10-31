pbdc_test_samtools = find_program('samtools')

pbdc_cram_tests = [
  'tiny',
]

test_env = [
  'ACTC=' + actc_main.full_path(),
  'MESON_BUILD_ROOT=' + meson.project_build_root(),
]

test(
  'pbdc formatting check',
  pbdc_clang_formatter,
  args : [
    '--all'],
  workdir : meson.project_source_root())

foreach t : pbdc_cram_tests
  test(
    'pbdc cram test - ' + t,
    pbdc_cram_script,
    args : [
      '--xunit-file=' + join_paths(meson.project_build_root(), 'pbdc-cram-' + t.underscorify() + '-test.xml'),
      '--verbose'] + files('cram/' + t + '.t'),
    env : test_env,
    is_parallel: not t.contains('pbdc'),
    timeout : 36000) # with '-O0 -g' tests can be *very* slow
endforeach
