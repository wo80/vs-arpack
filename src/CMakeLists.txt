include_directories(${f2c_INCLUDE_DIR})

# Add target: arpack-ng static library

add_library(arpack-ng STATIC)
add_library(ARPACK::arpack-ng ALIAS arpack-ng)

set_target_properties(arpack-ng PROPERTIES DEBUG_POSTFIX _d)

target_link_libraries(arpack-ng
  PUBLIC
  $<INSTALL_INTERFACE:$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,,LAPACK::LAPACK>>
  $<INSTALL_INTERFACE:$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,,BLAS::BLAS>>
  $<BUILD_INTERFACE:LAPACK::LAPACK>
  $<BUILD_INTERFACE:BLAS::BLAS>
)
#target_link_options(arpack-ng PUBLIC "${EXTRA_LDFLAGS}")
#set_target_properties(arpack-ng PROPERTIES VERSION 2.4.0)
#set_target_properties(arpack-ng PROPERTIES SOVERSION 2)
target_include_directories(arpack-ng
  PUBLIC
  # Exported location of headers
  $<INSTALL_INTERFACE:${ARPACK_INSTALL_INCLUDEDIR}>
  # Find arpackdef.h, arpackicb.h, stat*.h, debug*.h at build time
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  # For ICB interface
  #$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/ICB>
)

include_directories("${PROJECT_BINARY_DIR}") # Find arpackdef.h and arpackicb.h

add_subdirectory(arpack-ng)


# Add target: libarpack shared library

if(BUILD_SHARED_LIBS)
  add_library(libarpack SHARED)

  # Enable arcomplex type of ARPACK++
  target_compile_definitions(libarpack PRIVATE __GNUG__)

  target_link_libraries(libarpack arpack-ng ${SuperLU_LIBRARIES} BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})

  target_include_directories(libarpack
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
      ${CMAKE_CURRENT_SOURCE_DIR}/arpackpp/include
      ${CMAKE_CURRENT_SOURCE_DIR}/libarpack
      ${SuperLU_INCLUDE_DIR}
  )

  add_subdirectory(libarpack)
endif()
