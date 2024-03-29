#!/usr/bin/env bash
set -vex

#########
# BUILD #
#########

# configure
# '--wrap-mode nofallback' prevents meson from downloading
# stuff from the internet or using subprojects.
meson \
  --werror \
  --buildtype "${BUILDTYPE:-release}" \
  --default-library "${LIBRARYTYPE:-shared}" \
  --libdir lib \
  --unity "${ENABLED_UNITY_BUILD:-off}" \
  --prefix "${PREFIX_ARG:-/usr/local}" \
  -Db_coverage="${ENABLED_COVERAGE:-false}" \
  -Db_lto="${ENABLED_LTO:-false}" \
  -Db_sanitize="${ENABLED_SANITIZERS:-none}" \
  -Db_lundef="${ENABLED_LUNDEF:-true}" \
  -Dcpp_debugstl="${ENABLED_DEBUGSTL:-false}" \
  -Dtests="${ENABLED_TESTS:-false}" \
  -Dtests-internal="${ENABLED_INTERNAL_TESTS:-false}" \
  "${CURRENT_BUILD_DIR:-build}" .

#  --wrap-mode "${ENABLED_WRAP_MODE:-nofallback}" \
# build
ninja -C "${CURRENT_BUILD_DIR:-build}" -v
