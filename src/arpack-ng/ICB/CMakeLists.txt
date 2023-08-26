target_include_directories(arpack-ng
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_sources(arpack-ng
  PRIVATE
    arpack_icb.c
    debug.c
    stat.c)
