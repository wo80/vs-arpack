
target_include_directories(arpack-ng
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

add_subdirectory(SRC)
add_subdirectory(UTIL)
target_sources(arpack-ng PRIVATE f2c.c)

if(ENABLE_TESTS)
    add_subdirectory(EXAMPLES)
endif()
