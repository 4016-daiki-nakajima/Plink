add_executable(test_tetrahedralizer test/test_tetrahedralizer.cpp ${UTILS_SRC_FILES}
    ${GEOMETRY_SRCS}
)
target_link_libraries(test_tetrahedralizer PUBLIC igl::glfw igl_copyleft::tetgen)

add_executable(test_EigenSolver test/test_EigenSolver.cpp ${UTILS_SRC_FILES} ${GEOMETRY_SRCS})
target_link_libraries(test_EigenSolver PUBLIC igl::glfw igl_copyleft::tetgen)

add_executable(test_IO test/test_IO.cpp ${UTILS_SRC_FILES} ${GEOMETRY_SRCS})
target_link_libraries(test_IO PUBLIC igl::glfw igl_copyleft::tetgen)

add_executable(test_miniaudio test/test_miniaudio.cpp)
target_link_libraries(test_miniaudio PUBLIC miniaudio)