# 1_MyFirstCMakeProject

# Note: Everything after a # is a commentary!

cmake_minimum_required(VERSION 3.0) # CMake minimum version required. It is always good to ensure!

project (Project1) # name of your project

find_package(Eigen3 3.3 NO_MODULE) # tries to automatic find Eigen

#Catch the Eingen in the predefined path
if(NOT TARGET Eigen3::Eigen)
add_library(Eigen3::Eigen INTERFACE IMPORTED)
target_include_directories(Eigen3::Eigen INTERFACE /usr/local/include/eigen3)
endif()


include_directories(
	"${PROJECT_SOURCE_DIR}/FEM_GFEM"
    "${PROJECT_SOURCE_DIR}/FEM_GFEM/Inputs"
    "${PROJECT_SOURCE_DIR}/FEM_GFEM/Elements"
    "${PROJECT_SOURCE_DIR}/FEM_GFEM/lin_alg"
    "${PROJECT_SOURCE_DIR}/FEM_GFEM/Integration_points"
    "${PROJECT_SOURCE_DIR}/FEM_GFEM/Boundeary_conditions"
    "${PROJECT_SOURCE_DIR}/FEM_GFEM/shape_functions"
    "${PROJECT_SOURCE_DIR}/FEM_GFEM/generate_mesh"
    "${PROJECT_SOURCE_DIR}/FEM_GFEM/Post_processing"
)

add_subdirectory(FEM_GFEM)
add_subdirectory(third_party/external/matplotplusplus)

add_executable(myExecutable main.cpp) # creates an executable called myExecutable using the file main.cpp
target_link_libraries(myExecutable FEM_GFEM_library)
target_link_libraries(myExecutable matplot)