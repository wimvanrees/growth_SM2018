cmake_minimum_required(VERSION 3.14 FATAL_ERROR) 
project(ShellSolver LANGUAGES C CXX)

# location of additional cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# create json file with compilation output
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

## BASE DIRS
set(TBB_ROOT_DIR "/opt/homebrew/")
set(VTK_ROOT "/opt/homebrew/")
set(GSL_ROOT_DIR "/opt/homebrew/")

include(safeguards)
include(FetchContent)

# some definitions
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(TEST_SRC_DIR ${PROJECT_SOURCE_DIR}/test)
set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/external)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)


# fix in case different compiler is specified (http://stackoverflow.com/questions/29898219/why-does-cmake-ignore-addsystem-header-files-when-cxx-is-defined)
SET(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")


# compilation flags
message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
message(STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR}")
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS}  -fopenmp -Wall -Wcast-align -Wextra -Wformat-nonliteral -Wformat=2 -Winvalid-pch -Wmissing-format-attribute -Wmissing-include-dirs -Wswitch-default -Wswitch-enum -Wno-switch-default -Wno-missing-braces -Wno-switch-enum -Wno-suggest-attribute=format")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -march=native -mtune=native")
endif()


# c++14
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_EXTENSIONS off)

# find sources
file(GLOB SOURCES_SIM "${SRC_DIR}/simulations/*.cpp")
file(GLOB HEADERS_SIM "${SRC_DIR}/simulations/*.hpp")
file(GLOB SOURCES_LIB "${SRC_DIR}/libshell/*.cpp")
file(GLOB HEADERS_LIB "${SRC_DIR}/libshell/*.hpp")
file(GLOB SOURCES_HLBFGS "${EXT_PROJECTS_DIR}/hlbfgs/*.cpp")
file(GLOB SOURCES_TEST ${TEST_SRC_DIR}/testshell/*.cpp)

# ------------------------------------------------------------
# Create the main library
# ------------------------------------------------------------
add_library(libshell ${SOURCES_LIB} ${HEADERS_LIB})
target_include_directories(libshell PUBLIC "${SRC_DIR}/libshell/")

# ------------------------------------------------------------
# Main executable
# ------------------------------------------------------------
add_executable(shell ${SOURCES_SIM} ${HEADERS_SIM})
target_include_directories(shell PUBLIC "${SRC_DIR}/libshell/")
target_link_libraries(shell PUBLIC libshell)

# ------------------------------------------------------------
# DEPENDENCIES
# ------------------------------------------------------------

# OpenMP (uses the default FindOpenMP.cmake file from the CMake installation)
set(OpenMP_CXX_INCLUDE_DIR "${OPENMP_INCLUDE_PATH}")
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(libshell PUBLIC OpenMP::OpenMP_CXX)
else()
    message(WARNING 
        "Could not find an OpenMP library using CMake's FindOpenMP module. "
        "Not turning on OpenMP acceleration. "
    )
endif()

# EIGEN
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
message(STATUS "Found Eigen3 in: ${EIGEN3_INCLUDE_DIR}")
target_include_directories(libshell PUBLIC ${EIGEN3_INCLUDE_DIR})

# TRIANGLE (download directly from LIBIGL and create library)
FetchContent_Declare(
  triangle
  GIT_REPOSITORY https://github.com/libigl/triangle
  GIT_TAG master
)
FetchContent_MakeAvailable(triangle) # will try to find it locally first

message(STATUS "Downloaded Triangle in: ${triangle_SOURCE_DIR}")
target_compile_definitions(triangle PRIVATE TRILIBRARY)# -DNO_TIMER)

target_compile_definitions(libshell PUBLIC USETRILIBRARY)
target_include_directories(libshell PUBLIC ${triangle_SOURCE_DIR})
target_link_libraries(libshell PUBLIC triangle)

# HLBFGS (create library)
message(STATUS "HLBFGS")
add_library(hlbfgs ${SOURCES_HLBFGS})
target_include_directories(hlbfgs PRIVATE ${EIGEN3_INCLUDE_DIR})
target_compile_definitions(hlbfgs PRIVATE HLBFGS_EIGEN)

target_compile_definitions(libshell PUBLIC USEHLBFGS)
target_include_directories(libshell PUBLIC "${EXT_PROJECTS_DIR}/hlbfgs/")
target_link_libraries(libshell PUBLIC hlbfgs)

# LIBIGL
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG main)
FetchContent_MakeAvailable(libigl) # will try to find it locally first
message(STATUS "Downloaded IGL in: ${libigl_SOURCE_DIR}")
target_include_directories(libshell PUBLIC ${libigl_SOURCE_DIR}/include)

# TBB
set(TBB_USE_DEBUG_BUILD FALSE)
find_package(TBB REQUIRED)
target_include_directories(libshell PUBLIC ${TBB_INCLUDE_DIRS})
target_link_libraries(libshell PUBLIC TBB::tbb TBB::tbbmalloc)

# VTK (OPTIONAL BUT STRONGLY RECOMMENDED FOR IO)
find_package(VTK REQUIRED COMPONENTS CommonCore CommonDataModel RenderingCore FiltersCore FiltersSources IOXML)
if(VTK_FOUND)
    message(STATUS "Found VTK in: ${VTK_DIR}")
    target_compile_definitions(libshell PUBLIC USEVTK)
    if (VTK_MAJOR_VERSION VERSION_LESS "8.9.0")
        message(STATUS "This is an older version of VTK, using VTK_USE_FILE.")
        include(${VTK_USE_FILE})
        message(STATUS "Found vtk libraries: ${VTK_LIBRARIES}")
        target_link_libraries(libshell PUBLIC ${VTK_LIBRARIES})
    else()
        message(STATUS "This is a recent version of VTK, importing CMake targets from VTK build.")
    	target_link_libraries(libshell PUBLIC VTK::CommonCore VTK::CommonDataModel VTK::RenderingCore VTK::FiltersCore VTK::FiltersSources VTK::IOXML)
    endif()				   
else()
    message("Did not find VTK - compiling without")    
endif()

# GSL (OPTIONAL BUT SOME TESTS WONT WORK WITHOUT IT)
find_package(GSL)
if(GSL_FOUND)
    message(STATUS "Found GSL in: ${GSL_INCLUDE_DIRS}")
    target_include_directories(libshell PUBLIC ${GSL_INCLUDE_DIRS})    
    target_link_libraries(libshell PUBLIC ${GSL_LIBRARIES})
    target_compile_definitions(libshell PUBLIC USEGSL)
else()
    message(STATUS "Did not find GSL - compiling without")
endif()

# ------------------------------------------------------------
# Tests
# ------------------------------------------------------------
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG main
)
FetchContent_MakeAvailable(googletest)
SET(INSTALL_GTEST ON)
SET(BUILD_GMOCK OFF)

# testing
enable_testing()
add_executable(testshell ${SOURCES_TEST})
target_include_directories(testshell PUBLIC "${SRC_DIR}/libshell/")
target_link_libraries(testshell gtest gmock_main libshell)
add_test(NAME test1 COMMAND testshell)


