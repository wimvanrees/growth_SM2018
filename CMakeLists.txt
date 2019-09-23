CMAKE_MINIMUM_REQUIRED(VERSION 3.5)
set(PROJECT_NAME_STR ShellSolver)
project(${PROJECT_NAME_STR} C CXX)

# location of additional cmake modules
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    )

# create json file with compilation output
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(safeguards)

#https://github.com/Crascit/DownloadProject
include(cmake/DownloadProject.cmake)
if (CMAKE_VERSION VERSION_LESS 3.2)
    set(UPDATE_DISCONNECTED_IF_AVAILABLE "")
else()
    set(UPDATE_DISCONNECTED_IF_AVAILABLE "UPDATE_DISCONNECTED 1")
endif()

# some definitions
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(TEST_SRC_DIR ${PROJECT_SOURCE_DIR}/test)
set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/external)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)



# compilation flags
message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
message(STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR}")
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS}  -fopenmp -Wall -Wcast-align -Wextra -Wformat-nonliteral -Wformat=2 -Winvalid-pch -Wmissing-format-attribute -Wmissing-include-dirs -Wswitch-default -Wswitch-enum -Wno-switch-default -Wno-missing-braces -Wno-switch-enum -Wno-suggest-attribute=format")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -march=native -mtune=native")
endif()


# c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_EXTENSIONS off)

# find sources
file(GLOB SOURCES_SIM "${SRC_DIR}/simulations/*.cpp")
file(GLOB HEADERS_SIM "${SRC_DIR}/simulations/*.hpp")
file(GLOB SOURCES_LIB "${SRC_DIR}/libshell/*.cpp")
file(GLOB HEADERS_LIB "${SRC_DIR}/libshell/*.hpp")
file(GLOB SOURCES_HLBFGS "${EXT_PROJECTS_DIR}/hlbfgs/*.cpp")
file(GLOB SOURCES_TEST ${TEST_SRC_DIR}/testshell/*.cpp)

# compile main executable
add_executable(shell ${SOURCES_SIM} ${HEADERS_SIM})

# compile main library
add_library(libshell SHARED ${SOURCES_LIB} ${HEADERS_LIB})
target_link_libraries(shell libshell)

include_directories("${SRC_DIR}/libshell/")

# fix in case different compiler is specified (http://stackoverflow.com/questions/29898219/why-does-cmake-ignore-addsystem-header-files-when-cxx-is-defined)
SET(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")

# IGL
find_package(LIBIGL)
if(NOT LIBIGL_FOUND)
  download_project(PROJ               LIBIGL
                 GIT_REPOSITORY       https://github.com/libigl/libigl.git
                 GIT_TAG              master
                 ${UPDATE_DISCONNECTED_IF_AVAILABLE}
  )
  set(LIBIGL_INCLUDE_DIR "${LIBIGL_SOURCE_DIR}/include")
  message("Downloaded libIGL in: ${LIBIGL_SOURCE_DIR}")
else()
  message("Found libIGL in: ${LIBIGL_INCLUDE_DIR}")
endif()
include_directories(SYSTEM ${LIBIGL_INCLUDE_DIR}) # set to system to suppress deprecated warnings

# EIGEN
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
if(NOT EIGEN3_FOUND)
  set(EIGEN3_INCLUDE_DIR "${LIBIGL_SOURCE_DIR}/external/eigen/")	
  message("Using Eigen3 that comes with libIGL")
else()
  message("Found Eigen3 in: ${EIGEN3_INCLUDE_DIR}")
endif()
include_directories(SYSTEM "${EIGEN3_INCLUDE_DIR}") # set to system to suppress deprecated warnings

# TBB
# set variable to exclude debug versions
set(TBB_USE_DEBUG_BUILD FALSE)
find_package(TBB REQUIRED)
include_directories(${TBB_INCLUDE_DIRS})
target_link_libraries(libshell ${TBB_LIBRARIES})

# HLBFGS (create library)
include_directories("${EXT_PROJECTS_DIR}/HLBFGS/")
set(CMAKE_CXX_FLAGS    "${CMAKE_CXX_FLAGS} -DUSEHLBFGS -DHLBFGS_EIGEN")
# shared library
add_library(hlbfgs SHARED ${SOURCES_HLBFGS})
target_link_libraries(libshell hlbfgs)
target_link_libraries(shell hlbfgs)

# TRIANGLE (download directly from LIBIGL and create library)
download_project(PROJ		TRIANGLE
		GIT_REPOSITORY https://github.com/libigl/triangle
		GIT_TAG        d284c4a843efac043c310f5fa640b17cf7d96170
                 ${UPDATE_DISCONNECTED_IF_AVAILABLE}
)
include_directories(${TRIANGLE_SOURCE_DIR})
message("Included Triangle in: ${TRIANGLE_SOURCE_DIR}")	
set(CMAKE_CXX_FLAGS    "${CMAKE_CXX_FLAGS} -DUSETRILIBRARY")
add_library(triangle SHARED "${TRIANGLE_SOURCE_DIR}/triangle.c")
target_compile_definitions(triangle PRIVATE -DTRILIBRARY)# -DNO_TIMER)
target_link_libraries(shell triangle)
target_link_libraries(libshell triangle)

# VTK (OPTIONAL BUT STRONGLY RECOMMENDED FOR IO)
find_package(VTK)
if(VTK_FOUND)
    message("Found vtk in: ${VTK_INCLUDE_DIRS}")
    include(${VTK_USE_FILE})
    set(CMAKE_CXX_FLAGS    "${CMAKE_CXX_FLAGS} -DUSEVTK -DVTK_EXCLUDE_STRSTREAM_HEADERS")
    target_link_libraries(shell ${VTK_LIBRARIES})
    target_link_libraries(libshell ${VTK_LIBRARIES})
else()
    message("Did not find VTK - compiling without")    
endif()

# GSL (OPTIONAL BUT SOME TESTS WONT WORK WITHOUT IT)
find_package(GSL)
if(GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    set(CMAKE_CXX_FLAGS    "${CMAKE_CXX_FLAGS} -DUSEGSL")
    target_link_libraries(shell ${GSL_LIBRARIES})
    target_link_libraries(libshell ${GSL_LIBRARIES})
else()
    message("Did not find GSL - compiling without")
endif()

# DOXYGEN (https://tty1.net/blog/2014/cmake-doxygen_en.html)
find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install the HTML based API        
documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if(BUILD_DOCUMENTATION)
    if(NOT DOXYGEN_FOUND)
         message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif()

    set(doxyfile ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)

    message("Doxygen build started.")

    add_custom_target(doc
                      COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      COMMENT "Generating API documentation with Doxygen"
                      VERBATIM)

    #    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION     share/doc)
endif()


### FOR TESTING
download_project(PROJ                googletest
                 GIT_REPOSITORY      https://github.com/google/googletest.git
                 GIT_TAG             master
                 ${UPDATE_DISCONNECTED_IF_AVAILABLE}
)

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
# When using CMake 2.8.11 or later, header path dependencies
# are automatically added to the gtest and gmock targets.
# For earlier CMake versions, we have to explicitly add the
# required directories to the header search path ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include"
                        "${gmock_SOURCE_DIR}/include")
endif()

enable_testing()



# test example
set(PROJECT_TEST_NAME testshell)
add_executable(${PROJECT_TEST_NAME} ${SOURCES_TEST})
target_link_libraries(${PROJECT_TEST_NAME} gtest gmock_main libshell)
add_test(NAME test1 COMMAND ${PROJECT_TEST_NAME})



# make install
#include(install_target)
