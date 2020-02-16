cmake_minimum_required(VERSION 3.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# destination
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}")

# import boost library
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost REQUIRED)

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

include_directories(${Boost_INCLUDE_DIRS})

add_executable(gmsh2json gmsh2json.cpp)
target_link_libraries(gmsh2json ${Boost_LIBRARIES})
