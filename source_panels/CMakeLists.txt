cmake_minimum_required(VERSION 3.20)
project(source_panels)

# Add source files to the library
add_library(source_panels STATIC surface.cpp)

# Add third-party libraries
add_subdirectory(thirdparty/eigen)
add_subdirectory(thirdparty/rapidobj)

# Include directories for the library
target_include_directories(source_panels PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/eigen
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/rapidobj
)

# Link third-party libraries
target_link_libraries(source_panels PUBLIC Eigen3::Eigen rapidobj)




