#!/usr/bin/env bash
set -vex

export ENABLED_TESTS="true"
export ENABLED_INTERNAL_TESTS="${bamboo_ENABLED_INTERNAL_TESTS}"

case "${GCC_VERSION}" in
  next)
    module load gcc/8.1.0
    ;;

  clang)
    module load gtest/gcc48

    source /opt/rh/llvm-toolset-6.0/enable
    CC="clang"
    CXX="clang++"
    ;;

  *)
    module load gcc
    ;;
esac

module load ccache

export CC="ccache ${CC:-gcc}"
export CXX="ccache ${CXX:-g++}"
export CCACHE_BASEDIR="${PWD}"

# without -fno-sanitize-recover=all UBSAN failures won't abort the program
export CFLAGS="${CFLAGS} -fno-sanitize-recover=all"
export CXXFLAGS="${CXXFLAGS} -fno-sanitize-recover=all"

if [[ -z ${bamboo_planRepository_branchName+x} ]]; then
  : #pass
elif [[ ! -d /pbi/flash/bamboo/ccachedir ]]; then
  echo "[WARNING] /pbi/flash/bamboo/ccachedir is missing"
elif [[ $bamboo_planRepository_branchName == develop ]]; then
  export CCACHE_DIR=/pbi/flash/bamboo/ccachedir/${bamboo_shortPlanKey}.${bamboo_shortJobKey}.develop
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
elif [[ $bamboo_planRepository_branchName == master ]]; then
  export CCACHE_DIR=/pbi/flash/bamboo/ccachedir/${bamboo_shortPlanKey}.${bamboo_shortJobKey}.master
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
elif [[ $USER == bamboo ]]; then
  _shortPlanKey=$(echo ${bamboo_shortPlanKey}|sed -e 's/[0-9]*$//')
  export CCACHE_DIR=/pbi/flash/bamboo/ccachedir/${bamboo_shortPlanKey}.${bamboo_shortJobKey}
  if [[ -d /pbi/flash/bamboo/ccachedir/${_shortPlanKey}.${bamboo_shortJobKey}.develop ]]; then
    ( cd /pbi/flash/bamboo/ccachedir/
      cp -a ${_shortPlanKey}.${bamboo_shortJobKey}.develop $CCACHE_DIR
    )
  fi
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
fi
