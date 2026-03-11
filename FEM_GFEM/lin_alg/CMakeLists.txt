# 1_MyFirstCMakeProject

# Note: Everything after a # is a commentary!

#file(GLOB headers *.h) # creates a set headers with all the .h files
#file(GLOB sources *.cpp) # creates a set sources with all the .cpp files

# Create a set with what is inside headers and sources send it to parent folder 
#set(lin_alg_files ${headers} ${sources} PARENT_SCOPE) 

set(lin_alg_files
	${CMAKE_CURRENT_SOURCE_DIR}/lin_alg.h
	${CMAKE_CURRENT_SOURCE_DIR}/lin_alg.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Matrix.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Vector.cpp
	PARENT_SCOPE
)
