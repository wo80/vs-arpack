target_include_directories(arpack
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_sources(arpack
  PRIVATE
    arpack_icb.c
    debug.c
    stat.c)
