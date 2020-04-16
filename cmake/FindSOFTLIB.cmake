#.rst:
# FindSOFTLib
# -----------
#
# Finds SOFTLib, both binaries and headers.

set(SOFTLIB_INCLUDE_RELPATH1 "../../softlib/include/")
set(SOFTLIB_INCLUDE_RELPATH2 "../../softlib/build/include/")
set(SOFTLIB_LIBRARY_RELPATH "../../softlib/build/src/")

find_path(SOFTLIB_INCLUDE_DIR1 softlib/SFile.h
    HINTS ${SOFTLIB_INCLUDE_RELPATH1})
find_path(SOFTLIB_INCLUDE_DIR2 softlib/config.h
    HINTS ${SOFTLIB_INCLUDE_RELPATH2})
find_library(SOFTLIB_LIBRARY softlib libsoftlib
    HINTS ${SOFTLIB_LIBRARY_RELPATH})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SOFTLIB DEFAULT_MSG SOFTLIB_LIBRARY SOFTLIB_INCLUDE_DIR1 SOFTLIB_INCLUDE_DIR2)
mark_as_advanced(SOFTLIB_INCLUDE_DIR1 SOFTLIB_INCLUDE_DIR2 SOFTLIB_LIBRARY)

set(SOFTLIB_LIBRARIES ${SOFTLIB_LIBRARY})
set(SOFTLIB_INCLUDE_DIRS ${SOFTLIB_INCLUDE_DIR1} ${SOFTLIB_INCLUDE_DIR2})
