# masst_plus - Separate Mains
# Author: Yudong Liu <yudongl@andrew.cmu.edu>

# Set C++ standard
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Generate executable
add_executable(clustering_plus clustering_plus.cpp)
#target_sources(pairing_plus pairing_plus.cpp)

