cmake_minimum_required(VERSION 3.6)

project(ExaFMM_VPM)

if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(
    COMMAND git rev-parse --short HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH_SHORT
    OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
  set(GIT_BRANCH "")
  set(GIT_COMMIT_HASH "")
  set(GIT_COMMIT_HASH_SHORT "")
endif()

message(STATUS "On Git Branch: ${GIT_BRANCH} (${GIT_COMMIT_HASH})")

set(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(
    STATUS
      "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "${DEFAULT_BUILD_TYPE}"
      CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui, ccmake
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()
set(CMAKE_TRY_COMPILE_CONFIGURATION ${CMAKE_BUILD_TYPE})
set(CMAKE_CONFIGURATION_TYPES
    "Debug;Release"
    CACHE STRING "" FORCE)


enable_testing()

set(MSVC_COMPILER_FLAGS
    $<$<COMPILE_LANGUAGE:CXX>:/MP> # Use as many cores as available for
                                   # compilation.
    $<$<COMPILE_LANGUAGE:CXX>:/W4> # /Wall # Baseline reasonable warnings.
    $<$<COMPILE_LANGUAGE:CXX>:/permissive-> # Standards conformance mode for
                                            # MSVC compiler.
    $<$<COMPILE_LANGUAGE:CXX>:$<$<CONFIG:Release>:/Oi>>
    $<$<COMPILE_LANGUAGE:CXX>:$<$<CONFIG:Release>:/O2>>
    $<$<COMPILE_LANGUAGE:CXX>:$<$<CONFIG:Release>:/Ot>>
    $<$<COMPILE_LANGUAGE:CXX>:$<$<CONFIG:Release>:/GL>>
    $<$<COMPILE_LANGUAGE:CXX>:$<$<CONFIG:Release>:/GF>>)

set(CLANG_COMPILER_FLAGS
    -g
    -W
    # -Wall
    -Wextra # Reasonable and standard.
    -Wshadow # Warn the user if a variable declaration shadows one from a parent
             # context.
    -Wnon-virtual-dtor # Warn the user if a class with virtual functions has a
                       # non-virtual destructor. This helps catch hard to track
                       # down memory errors.
    -Wcast-align # Warn for potential performance problem casts.
    -Wunused # Warn on anything being unused.
    -Wunused-function
    -Wunused-variable
    -Wunused-parameter
    -Woverloaded-virtual # Warn if you overload (not override) a virtual
                         # function.
    $<$<CONFIG:Debug>:-O0;-ggdb;-fno-elide-constructors>)

set(GCC_COMPILER_FLAGS
    ${CLANG_COMPILER_FLAGS}
)

if(MSVC)
  set(COMPILER_FLAGS ${MSVC_COMPILER_FLAGS})
elseif(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
  set(COMPILER_FLAGS ${CLANG_COMPILER_FLAGS})
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(COMPILER_FLAGS ${GCC_COMPILER_FLAGS})
else()
  message(
    AUTHOR_WARNING
      "No compiler flags set for '${CMAKE_CXX_COMPILER_ID}' compiler.")
endif()



find_package(OpenMP)

# Add executables for exafmm
add_executable(exafmm exafmm_vpm/fmm.cpp)
add_executable(exafmm_kernel exafmm_vpm/kernel.cpp)

# Include directories for exafmm
target_include_directories(exafmm PRIVATE exafmm_vpm)
target_include_directories(exafmm_kernel PRIVATE exafmm_vpm)

# Add source_panels as a subdirectory
add_subdirectory(source_panels)

# Link the source_panels library to exafmm and exafmm_kernel using keyword signature
target_link_libraries(exafmm PRIVATE source_panels)
target_link_libraries(exafmm_kernel PRIVATE source_panels)

# OpenMP settings
if (OpenMP_CXX_FOUND)
  target_link_libraries(exafmm PRIVATE OpenMP::OpenMP_CXX)
  target_link_libraries(exafmm_kernel PRIVATE OpenMP::OpenMP_CXX)
endif()

# Compile options
target_compile_options(exafmm PRIVATE ${COMPILER_FLAGS})
target_compile_options(exafmm_kernel PRIVATE ${COMPILER_FLAGS})

# Compile definitions
target_compile_definitions(exafmm PRIVATE EXAFMM_EAGER)

