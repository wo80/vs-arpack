include_directories(${f2c_INCLUDE_DIR})

# Add target: arpack static library

add_library(arpack STATIC)
add_library(ARPACK::arpack ALIAS arpack)

set_target_properties(arpack PROPERTIES DEBUG_POSTFIX _d)

target_link_libraries(arpack
  PUBLIC
  $<INSTALL_INTERFACE:$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,,LAPACK::LAPACK>>
  $<INSTALL_INTERFACE:$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,,BLAS::BLAS>>
  $<BUILD_INTERFACE:LAPACK::LAPACK>
  $<BUILD_INTERFACE:BLAS::BLAS>
)
#target_link_options(arpack PUBLIC "${EXTRA_LDFLAGS}")
#set_target_properties(arpack PROPERTIES VERSION 2.4.0)
#set_target_properties(arpack PROPERTIES SOVERSION 2)
target_include_directories(arpack
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


# Add target: arpack-ds shared library

add_library(arpack-ds SHARED)

if(BUILD_SHARED_LIBS)
  add_library(ARPACK::arpack-ds ALIAS arpack-ds)

  # Enable arcomplex type of ARPACK++
  target_compile_definitions(arpack-ds PRIVATE __GNUG__)

  target_link_libraries(arpack-ds arpack superlu::superlu BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})

  target_include_directories(arpack-ds
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/arpackpp/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/arpack-ds>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

  add_subdirectory(arpack-ds)
endif()
