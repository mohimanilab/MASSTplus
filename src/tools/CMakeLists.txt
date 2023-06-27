# masst_plus - Main
# Author: Yudong Liu <yudongl@andrew.cmu.edu>

# Add subdirectories
add_subdirectory(masst_plus)
add_subdirectory(pairing_plus)
add_subdirectory(clustering_plus)
add_subdirectory(networking_plus)

# Link libraries
target_link_libraries(masst_plus 
	MGF)