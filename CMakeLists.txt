cmake_minimum_required(VERSION 3.2.0)
include(CheckCXXCompilerFlag)

set(GT_FOLDER "${CMAKE_SOURCE_DIR}/../GeometricTools/GTE")

project(aabbcyl)

set (CMAKE_VERBOSE_MAKEFILE ON)

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Debug|Release|RelWithDebInfo|MinSizeRel" FORCE)
endif()

check_cxx_compiler_flag (-Wall HAS_WALL)
if (HAS_WALL)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

#find_package (OpenMP)
#if (OpenMP_CXX_FOUND)
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif()

check_cxx_compiler_flag (-xHost HAS_XHOST)
check_cxx_compiler_flag (-march=native HAS_MARCH_NATIVE)
if (HAS_XHOST)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost")
elseif (HAS_MARCH_NATIVE)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

check_cxx_compiler_flag (-ipo HAS_IPO)
if (HAS_IPO)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ipo")
endif()

check_cxx_compiler_flag (-flto HAS_FLTO)
if (HAS_FLTO)
    set (CMAKE_CXX_FLAGS "-flto ${CMAKE_CXX_FLAGS}")
endif()

enable_testing()

add_executable(test1 test1.cpp)
add_executable(test2 test2.cpp)
add_executable(test3 test3.cpp)

target_include_directories(test1 PRIVATE ${GT_FOLDER})
target_include_directories(test2 PRIVATE ${GT_FOLDER})
target_include_directories(test3 PRIVATE ${GT_FOLDER})

add_test(test1 test1)
add_test(test2 test2)
add_test(test3 test3)


