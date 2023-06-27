# masst_plus

cmake_minimum_required (VERSION 2.8)

# set the project name
project(masst_plus)

# set C++ standard

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-static-libstdc++")

# filesystem support
link_libraries(stdc++fs)

# external dependencies
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external)

# Add subdirectories
add_subdirectory(src)
