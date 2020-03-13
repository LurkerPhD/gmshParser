cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# destination
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

# gmsh {{{
find_library(GMSH_LIB gmsh)
message("there is the gmsh lib :" ${GMSH_LIB})
if(NOT GMSH_LIB)
  message(FATAL_ERROR "Could not find libgmsh")
endif(NOT GMSH_LIB)

find_path(GMSH_INC gmsh.h)
message("there is the gmsh.h :" ${GMSH_INC})
if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find gmsh.h")
endif(NOT GMSH_INC)

if(GMSH_LIB MATCHES ".a") # FIXME - generalize this
  find_library(BLAS_LIB blas)
  if(BLAS_LIB)
    list(APPEND EXTRA_LIBS ${BLAS_LIB})
  endif(BLAS_LIB)
  find_library(LAPACK_LIB lapack)
  if(LAPACK_LIB)
    list(APPEND EXTRA_LIBS ${LAPACK_LIB})
  endif(LAPACK_LIB)
endif(GMSH_LIB MATCHES ".a")

if(WIN32 OR CYGWIN)
  list(APPEND EXTRA_LIBS winmm wsock32 ws2_32 psapi)
endif(WIN32 OR CYGWIN)
# }}}

# import boost library
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost REQUIRED)

include_directories(${GMSH_INC})
include_directories(${Boost_INCLUDE_DIRS})

add_executable(gmshParser gmsh2json.cpp)
target_link_libraries(gmshParser ${Boost_LIBRARIES})
target_link_libraries(gmshParser ${GMSH_LIB})
