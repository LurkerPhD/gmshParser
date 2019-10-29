cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# TODO: Gmsh should create a cmake file when installing the lib, so all the
# external dependencies can be obtained automatically

set(CMAKE_CXX_STANDARD 17)

project(gmshParser)

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

include_directories(${GMSH_INC})

# add_executable(${PROJECT_NAME} gmshParser.cpp tinyxml2.cpp)
add_executable(${PROJECT_NAME} gmshParser.cpp main.cpp tinyxml2.cpp)
target_link_libraries(${PROJECT_NAME} ${GMSH_LIB} ${EXTRA_LIBS})

set(POST_PROCESSER gmshPostProcess)
add_executable(${POST_PROCESSER} gmshPostProcess.cpp)
target_link_libraries(${POST_PROCESSER} ${GMSH_LIB} ${EXTRA_LIBS})
