project(Vaango)
cmake_minimum_required(VERSION 2.8)

# Add requirements for PCL
find_package(PCL 1.2 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

# Create a library for the PLY converter
add_executable(PLYFileReader PLYFileReader.cc)
target_link_libraries (PLYFileReader ${PCL_LIBRARIES})
