# Copyright (c) 2026 Huawei Technologies Co., Ltd. All Rights Reserved.
#
# KuQCD is licensed under Mulan PSL v2.
# You can use this software according to the terms and conditions of the Mulan PSL v2.
# You may obtain a copy of Mulan PSL v2 at:
#       http://license.coscl.org.cn/MulanPSL2
# THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT, MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE.
# See the Mulan PSL v2 for more details.

####################################################################################
# basic set up
####################################################################################
cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# Disable gnu exentions
set(CMAKE_CXX_EXTENSIONS OFF)

# Disable in source builds this is only a temporary fix, but for now we need it as
# cmake will otherwise overwrite the existing makefiles.
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
# add a directory for cmake modules

project("KuQCD" VERSION 0.0)

# Print the configuration details to stdout
message(STATUS "")
message(STATUS "${PROJECT_NAME} ${PROJECT_VERSION} **")
message(STATUS "cmake version: ${CMAKE_VERSION}")
message(STATUS "Source location: ${CMAKE_SOURCE_DIR}")
message(STATUS "Build location: ${CMAKE_BINARY_DIR}")

include_directories(${CMAKE_SOURCE_DIR}/include)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)


# ######################################################################################################################
# basic options
# ######################################################################################################################
# option(USE_GCC "Use GCC compiler" ON)
# option(USE_CLANG "Use clang compiler" OF)

set(KUNPENG_ARCH "SVE512" CACHE STRING "set kunpeng architecture (SVE256 or SVE512)")
set_property(CACHE KUNPENG_ARCH PROPERTY STRINGS "SVE256" "SVE512")

option(DEBUG "Debug" OFF)
option(USE_SVE "Use arm SVE" ON)
option(USE_HUGE_PAGE_MEMORY "Use huge page memory" OFF)
option(USE_KUPL "Use KUPL" ON)
option(USE_KUPL_COMM "Use KUPL for communication" ON)
option(USE_KUPL_ALL "Use KUPL for bulk data" OFF)
option(USE_KPMU "Use KPMU" OFF)
option(USE_JEMALLOC "Use jemalloc" OFF)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_C_COMPILER clang)
    set(CMAKE_CXX_COMPILER clang++)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_C_COMPILER gcc)
    set(CMAKE_CXX_COMPILER g++)
else()
    message(FATAL_ERROR "invalid COMPLIER")
endif()
set(CMAKE_SKIP_RPATH TRUE)

if(DEBUG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fstack-protector-strong")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fstack-protector-strong")

set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -s -Wl,-z,relro -Wl,-z,now -Wl,-z,noexecstack")

if(USE_KUPL)
    if(USE_KUPL_COMM AND USE_KUPL_ALL)
        message(SEND_ERROR "cannot use KUPL shared memory for both surface and bulk")
    endif()
    add_definitions(-DUSE_KUPL)
    if(USE_KUPL_COMM)
        add_definitions(-DUSE_KUPL_COMM)
    endif()
    if(USE_KUPL_ALL)
        add_definitions(-DUSE_KUPL_ALL)
    endif()
endif()

if(USE_HUGE_PAGE_MEMORY)
    add_definitions(-DUSE_HPMEM)
endif()


# ######################################################################################################################
# advanced options
# ######################################################################################################################
option(KUQCD_BUILD_ALL_TESTS "build tests by default" ON)
option(KUQCD_INSTALL_ALL_TESTS "install tests by default" ON)
option(KUQCD_BUILD_SHAREDLIB "build kuqcd as a shared lib" ON)
option(KUQCD_MPI "enable MPI" ON)
option(KUQCD_OPENMP "enable OpenMP" ON)

set(KUQCD_CXX_STANDARD 17 CACHE STRING "set the CXX Standard (11, 14 or 17)")
set_property(CACHE KUQCD_CXX_STANDARD PROPERTY STRINGS 11 14 17)

mark_as_advanced(KUQCD_BUILD_ALL_TESTS)
mark_as_advanced(KUQCD_INSTALL_ALL_TESTS)

####################################################################################
# configuration set up
####################################################################################

### bisheng5.0 need
add_compile_options(-Werror -Wall -Wno-instantiation-after-specialization -fPIC -fstack-protector -fvisibility=hidden -D_FORTIFY_SOURCE=2 -ftrapv)
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 19.0)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wno-missing-designated-field-initializers>)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wno-vla-cxx-extension>)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-mllvm>)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-treat-scalable-fixed-error-as-warning=true>)
endif()
else()
add_compile_options(-Wno-array-bounds -Wno-maybe-uninitialized -Wno-unknown-pragmas)
endif()

set(CPU_ARCH ${CMAKE_HOST_SYSTEM_PROCESSOR})
if(${CPU_ARCH} STREQUAL "aarch64")
    if(${KUNPENG_ARCH} STREQUAL "SVE512")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mcpu=hip11 -ffast-math -funroll-loops -ftree-vectorize -fopenmp")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=hip11 -ffast-math -funroll-loops -ftree-vectorize -fopenmp")
        if(USE_SVE)
            set(SVEBITS 512)
            message("use sve 512")
        endif()
    elseif(${KUNPENG_ARCH} STREQUAL "SVE256")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mcpu=hip09 -ffast-math -funroll-loops -ftree-vectorize -fopenmp")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=hip09 -ffast-math -funroll-loops -ftree-vectorize -fopenmp")
        if(USE_SVE)
            set(SVEBITS 256)
            message("use sve 256")
        endif()
    else()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mcpu=native -ffast-math -funroll-loops -ftree-vectorize -fopenmp")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=native -ffast-math -funroll-loops -ftree-vectorize -fopenmp")
    endif()
elseif(${CPU_ARCH} STREQUAL "x86_64")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=native -funroll-loops -ftree-vectorize -fopenmp")
elseif(${CPU_ARCH} STREQUAL "ppc64le")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -mcpu=native -mtune=native -funroll-loops -ftree-vectorize -fopenmp")
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

message(STATUS "${CMAKE_CXX_FLAGS}")

set(CMAKE_CXX_STANDARD ${KUQCD_CXX_STANDARD})

add_subdirectory(src)
# add_subdirectory(thirdparty)

# file(GLOB_RECURSE RESOURCES "${CMAKE_SOURCE_DIR}/testSrc/*")
