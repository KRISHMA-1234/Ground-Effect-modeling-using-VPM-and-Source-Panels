# Install script for directory: /home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/AdolcForward"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/AlignedVector3"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/ArpackSupport"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/AutoDiff"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/BVH"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/EulerAngles"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/FFT"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/IterativeSolvers"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/KroneckerProduct"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/LevenbergMarquardt"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/MatrixFunctions"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/MPRealSupport"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/NNLS"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/NonLinearOptimization"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/NumericalDiff"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/OpenGLSupport"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/Polynomials"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/SparseExtra"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/SpecialFunctions"
    "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/Splines"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "/home/krishma/blade_particle/fmm_code/to_karishma/source_panels/thirdparty/eigen/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/krishma/blade_particle/fmm_code/to_karishma/build/source_panels/thirdparty/eigen/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

