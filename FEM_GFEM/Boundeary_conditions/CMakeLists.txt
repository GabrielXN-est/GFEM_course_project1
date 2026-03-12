# 1_MyFirstCMakeProject

# Note: Everything after a # is a commentary!

#file(GLOB headers *.h) # creates a set headers with all the .h files
#file(GLOB sources *.cpp) # creates a set sources with all the .cpp files

# Create a set with what is inside headers and sources send it to parent folder 
#set(BC ${headers} ${sources} PARENT_SCOPE) 

set(BC
	${CMAKE_CURRENT_SOURCE_DIR}/body_func.h
	${CMAKE_CURRENT_SOURCE_DIR}/Bondeary_conditions.h
	PARENT_SCOPE
)