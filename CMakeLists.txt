# 1_MyFirstCMakeProject

# Note: Everything after a # is a commentary!

cmake_minimum_required(VERSION 3.0) # CMake minimum version required. It is always good to ensure!

project (Project1) # name of your project

include_directories(
	"${PROJECT_SOURCE_DIR}/FEM_GFEM"
    "${PROJECT_SOURCE_DIR}/FEM_GFEM/Inputs"
    "${PROJECT_SOURCE_DIR}/FEM_GFEM/Elements"
    "${PROJECT_SOURCE_DIR}/FEM_GFEM/lin_alg"
    "${PROJECT_SOURCE_DIR}/FEM_GFEM/Integration_points"
    "${PROJECT_SOURCE_DIR}/FEM_GFEM/Boundeary_conditions"
    "${PROJECT_SOURCE_DIR}/FEM_GFEM/shape_functions"
)

add_subdirectory(FEM_GFEM)

add_executable(myExecutable main.cpp) # creates an executable called myExecutable using the file main.cpp
target_link_libraries(myExecutable library)