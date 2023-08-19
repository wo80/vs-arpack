include(FindPackageHandleStandardArgs)

if (NOT WIN32)
    find_package(PkgConfig)
    if (PKG_CONFIG_FOUND)
         pkg_check_modules(PKG_LIBF2C f2c)
    endif ()
endif (NOT WIN32)

find_path(f2c_INCLUDE_DIR f2c.h
    ${PKG_LIBF2C_INCLUDE_DIRS}
    /usr/include
    /usr/local/include
)

find_library(f2c_LIBRARY
    NAMES
    f2c libf2c
    PATHS
    ${PKG_LIBF2C_LIBRARY_DIRS}
    /usr/lib
    /usr/local/lib
)

find_package_handle_standard_args(F2C DEFAULT_MSG f2c_LIBRARY f2c_INCLUDE_DIR)
