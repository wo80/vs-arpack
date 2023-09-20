
target_include_directories(arpack
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

add_subdirectory(SRC)
add_subdirectory(UTIL)

if(ENABLE_ICB)
  add_subdirectory(ICB)
endif()

target_sources(arpack PRIVATE f2c.c)

if(ENABLE_TESTS)
  add_subdirectory(EXAMPLES)
endif()
