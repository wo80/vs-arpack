include(FindPackageHandleStandardArgs)

if (NOT WIN32)
    find_package(PkgConfig)
    if (PKG_CONFIG_FOUND)
         pkg_check_modules(PKG_LIBSUPERLU libsuperlu)
    endif ()
endif (NOT WIN32)

find_path(SuperLU_INCLUDE_DIR SuperLU/slu_util.h
    ${PKG_LIBSUPERLU_INCLUDE_DIRS}
    /usr/include
    /usr/local/include
)

find_library(SuperLU_LIBRARIES
    NAMES
    superlu superlu.dll
    PATHS
    ${PKG_LIBSUPERLU_LIBRARY_DIRS}
    /usr/lib
    /usr/local/lib
)

find_package_handle_standard_args(SuperLU DEFAULT_MSG SuperLU_LIBRARIES SuperLU_INCLUDE_DIR)
