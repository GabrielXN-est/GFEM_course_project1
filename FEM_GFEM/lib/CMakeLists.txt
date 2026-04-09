# 1_MyFirstCMakeProject

add_library(FEM_GFEM_library ${base_files} ${element_files} ${lin_alg_files} ${input_files} ${integration_points_files} ${BC} ${shape_functions} ${post_process})
target_link_libraries(FEM_GFEM_library PUBLIC Eigen3::Eigen)