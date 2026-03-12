# 1_MyFirstCMakeProject

# Note: Everything after a # is a commentary!

#file(GLOB headers *.h) # creates a set headers with all the .h files
#file(GLOB sources *.cpp) # creates a set sources with all the .cpp files

# Create a set with what is inside headers and sources send it to parent folder 
#set(element_files ${headers} ${sources} PARENT_SCOPE) 

set(element_files
	${CMAKE_CURRENT_SOURCE_DIR}/bars.h
	${CMAKE_CURRENT_SOURCE_DIR}/bars.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/node.h
	${CMAKE_CURRENT_SOURCE_DIR}/element.h
	${CMAKE_CURRENT_SOURCE_DIR}/element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bar.cpp
	PARENT_SCOPE
)