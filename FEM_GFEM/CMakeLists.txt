# 1_MyFirstCMakeProject

# Note: Everything after a # is a commentary!

file(GLOB headers *.h) # creates a set headers with all the .h files
file(GLOB sources *.cpp) # creates a set sources with all the .cpp files

# Create a set with what is inside headers and sources send it to parent folder 
set(base_files ${headers} ${sources} PARENT_SCOPE) 

add_subdirectory(Inputs) 
add_subdirectory(Elements)
add_subdirectory(lin_alg)
add_subdirectory(Integration_points)
add_subdirectory(lib)