# masst_plus - Separate Mains
# Author: Yudong Liu <yudongl@andrew.cmu.edu>

# Generate executable
add_executable(load load_library.cpp ../../spectral_index.hpp ../spectral_index.cpp)
add_executable(search search.cpp ../../spectral_index.hpp ../spectral_index.cpp)

# Link libraries
target_link_libraries(load
	MGF)
target_link_libraries(search
	MGF)
