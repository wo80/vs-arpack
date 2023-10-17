cmake_minimum_required(VERSION 3.16.0)

project(arpack VERSION 3.9.1
               DESCRIPTION "ARPACK-NG: solve large scale eigenvalue problems"
               LANGUAGES C CXX)

set(arpack_ng_MAJOR_VERSION 3)
set(arpack_ng_MINOR_VERSION 9)
set(arpack_ng_PATCH_VERSION 1)
set(arpack_ng_VERSION ${arpack_ng_MAJOR_VERSION}.${arpack_ng_MINOR_VERSION}.${arpack_ng_PATCH_VERSION})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# set C99 standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(BUILD_SHARED_LIBS "Build shared libraries instead of static libraries" ON)
option(USE_EIGEN "Enable support for eigenvalue-problems solver based on ICB and eigen" OFF)
option(ENABLE_TESTS "Compile ARPACK tests" ON)
option(ENABLE_ICB "Include ICB interface" OFF)
option(ENABLE_TIMING "Enable ARPACK timing feature" ON)
option(ENABLE_TRACING "Enable ARPACK tracing feature" ON)

# 'make install' to the correct location, and also define
# paths for target_include_directories and target_link_libraries
include(GNUInstallDirs)
set(ARPACK_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/arpack")
set(ARPACK_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/arpackng")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

if(MSVC)
  add_compile_definitions(NOMINMAX _COMPLEX_DEFINED)
endif()

if(NOT ENABLE_TIMING)
  add_compile_definitions(NO_TIMER)
endif()

if(NOT ENABLE_TRACING)
  add_compile_definitions(NO_TRACE)
endif()

if(BUILD_SHARED_LIBS)
  find_package(superlu REQUIRED)

  if (TARGET superlu::superlu)
    message( STATUS "Found SuperLU: version ${superlu_VERSION}" )
  else()
    message(FATAL_ERROR "SuperLU not found")
  endif()
endif()

# Find BLAS

if (NOT TARGET BLAS::BLAS) # Search only if not already found by upper CMakeLists.txt
  find_package(BLAS REQUIRED)

  # BLAS::BLAS target was already created at this point by FindBLAS.cmake if cmake version >= 3.18
  if (NOT TARGET BLAS::BLAS) # Create target "at hand" to ensure compatibility if cmake version < 3.18
    add_library(BLAS::BLAS INTERFACE IMPORTED)
    set_target_properties(BLAS::BLAS PROPERTIES INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}")
  endif()
endif()

# Find LAPACK

if (NOT TARGET LAPACK::LAPACK) # Search only if not already found by upper CMakeLists.txt
  find_package(LAPACK REQUIRED)

  # LAPACK::LAPACK target was already created at this point by FindLAPACK.cmake if cmake version >= 3.18
  if (NOT TARGET LAPACK::LAPACK) # Create target "at hand" to ensure compatibility if cmake version < 3.18
    add_library(LAPACK::LAPACK INTERFACE IMPORTED)
    set_target_properties(LAPACK::LAPACK PROPERTIES INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
  endif()
endif()

# Find Eigen3

if (USE_EIGEN)
  find_package(Eigen3 3.3 QUIET)
  if (NOT EIGEN3_FOUND) # If not found, piggy-back pkg-config files.
    message(WARNING "CMake didn't find the Eigen3 package. Try to look for pkg-config file...")
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(EIGEN3 REQUIRED eigen3>=3.3)
    set(EIGEN3_INCLUDE_DIR ${EIGEN3_INCLUDE_DIRS})
  endif ()
  find_program (BASH_PROGRAM bash)

  # Look for headers.
  find_path(EIGEN3_ITERATIVE_SOLVER_DIR NAMES IterativeLinearSolvers PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
  if (NOT EIGEN3_ITERATIVE_SOLVER_DIR)
    message(FATAL_ERROR "-- Eigen/IterativeLinearSolvers header not found.")
  endif ()
  find_path(EIGEN3_SPARSE_DIR NAMES Sparse PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
  if (NOT EIGEN3_SPARSE_DIR)
    message(FATAL_ERROR "-- Eigen/Sparse header not found.")
  endif ()
  find_path(EIGEN3_SPARSELU_SOLVER_DIR NAMES SparseLU PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
  if (NOT EIGEN3_SPARSELU_SOLVER_DIR)
    message(FATAL_ERROR "-- Eigen/SparseLU header not found.")
  endif ()
  find_path(EIGEN3_SPARSEQR_SOLVER_DIR NAMES SparseQR PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
  if (NOT EIGEN3_SPARSEQR_SOLVER_DIR)
    message(FATAL_ERROR "-- Eigen/SparseQR header not found.")
  endif ()
  find_path(EIGEN3_SPARSECHOLESKY_SOLVER_DIR NAMES SparseCholesky PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
  if (NOT EIGEN3_SPARSECHOLESKY_SOLVER_DIR)
    message(FATAL_ERROR "-- Eigen/SparseCholesky header not found.")
  endif ()
  find_path(EIGEN3_DENSE_DIR NAMES Dense PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
  if (NOT EIGEN3_DENSE_DIR)
    message(FATAL_ERROR "-- Eigen/Dense header not found.")
  endif ()
  find_path(EIGEN3_DENSELU_DIR NAMES LU PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
  if (NOT EIGEN3_DENSELU_DIR)
    message(FATAL_ERROR "-- Eigen/LU header not found.")
  endif ()
  find_path(EIGEN3_DENSEQR_DIR NAMES QR PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
  if (NOT EIGEN3_DENSEQR_DIR)
    message(FATAL_ERROR "-- Eigen/QR header not found.")
  endif ()
  find_path(EIGEN3_DENSECHOLESKY_DIR NAMES Cholesky PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
  if (NOT EIGEN3_DENSECHOLESKY_DIR)
    message(FATAL_ERROR "-- Eigen/Cholesky header not found.")
  endif ()

  set(ICB "ON")
endif ()

if(ENABLE_TESTS)
  enable_testing()
  # set(CMAKE_CTEST_COMMAND ctest -V)
endif()

add_subdirectory(src)


# Install

if(BUILD_SHARED_LIBS)
  install(TARGETS arpack arpack-ds
    EXPORT arpackngTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

install(FILES
  src/arpack-ng/arpack.h
  src/arpack-ng/arpack_types.h
DESTINATION
  ${ARPACK_INSTALL_INCLUDEDIR})
    
if(ICB)
  install(FILES
    src/arpack-ng/ICB/arpack.hpp
    #src/arpack-ng/ICB/debug_c.h
    #src/arpack-ng/ICB/debug_c.hpp
    #src/arpack-ng/ICB/stat_c.h
    #src/arpack-ng/ICB/stat_c.hpp
  DESTINATION
    ${ARPACK_INSTALL_INCLUDEDIR})

  #if (EIGEN)
  #  install(FILES EXAMPLES/MATRIX_MARKET/arpackSolver.hpp DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
  #endif()
endif()

#configure_file(arpackdef.h.in "${PROJECT_BINARY_DIR}/arpackdef.h" @ONLY)
#install(FILES "${PROJECT_BINARY_DIR}/arpackdef.h" DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})

install(EXPORT arpackngTargets
  DESTINATION "${ARPACK_INSTALL_CMAKEDIR}"
)
# Provide find_package for arpack-ng to users.
configure_file(cmake/arpackng-config.cmake.in "${PROJECT_BINARY_DIR}/arpackng-config.cmake" @ONLY)
configure_file(cmake/arpackng-config-version.cmake.in "${PROJECT_BINARY_DIR}/arpackng-config-version.cmake" @ONLY)
install(
  FILES
  "${PROJECT_BINARY_DIR}/arpackng-config.cmake"
  "${PROJECT_BINARY_DIR}/arpackng-config-version.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/arpackng) # find_package(arpackng)
