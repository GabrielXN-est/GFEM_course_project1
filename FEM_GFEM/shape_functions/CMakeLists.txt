# 1_MyFirstCMakeProject

# Note: Everything after a # is a commentary!

#file(GLOB headers *.h) # creates a set headers with all the .h files
#file(GLOB sources *.cpp) # creates a set sources with all the .cpp files

# Create a set with what is inside headers and sources send it to parent folder 
#set(base_files ${headers} ${sources} PARENT_SCOPE) 

# Base FEM sources used by the library
set(shape_functions
	${CMAKE_CURRENT_SOURCE_DIR}/shape_functions.h
	${CMAKE_CURRENT_SOURCE_DIR}/shape_functions.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Enrichment.h
	${CMAKE_CURRENT_SOURCE_DIR}/Enrichment.cpp
	PARENT_SCOPE
)