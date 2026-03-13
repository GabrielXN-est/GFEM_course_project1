# 1_MyFirstCMakeProject

# Note: Everything after a # is a commentary!

#file(GLOB headers *.h) # creates a set headers with all the .h files
#file(GLOB sources *.cpp) # creates a set sources with all the .cpp files

# Create a set with what is inside headers and sources send it to parent folder 
#set(base_files ${headers} ${sources} PARENT_SCOPE) 

# Base FEM sources used by the library
set(Generate_mesh_files
	${CMAKE_CURRENT_SOURCE_DIR}/create_input.h
	${CMAKE_CURRENT_SOURCE_DIR}/create_input.cpp
)

add_library(mesh_library ${Generate_mesh_files})