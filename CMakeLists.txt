cmake_minimum_required(VERSION 3.10)
project(SimpleModal)

set(CMAKE_CXX_STANDARD 17)
# cmake -B build -S .
# -DCMAKE_TOOLCHAIN_FILE=C:/Users/Zhehao/Documents/Libraries/vcpkg/scripts/buildsystems/vcpkg.cmake
# set(CMAKE_TOOLCHAIN_FILE
# "C:/Users/Zhehao/Documents/Libraries/vcpkg/scripts/buildsystems/vcpkg.cmake")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Define project root directory
set(PROJECT_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
add_definitions(-DPROJECT_ROOT_DIR="${PROJECT_ROOT_DIR}")

set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS TRUE)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# find_package(fmt CONFIG REQUIRED) find_package(glfw3 3.4 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(Eigen3 REQUIRED)
# find_package(Boost REQUIRED)

message(CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE})
# find_package(PkgConfig REQUIRED) pkg_check_modules(gmp REQUIRED
# IMPORTED_TARGET gmp) find_package(GMP REQUIRED)

# add_subdirectory(external/fTetWild)

add_subdirectory(external/libigl)
igl_include(glfw)
igl_include(opengl)
igl_include(imgui)
igl_include(embree)
igl_include(copyleft tetgen)

add_subdirectory(external/spectra)
add_subdirectory(external/miniaudio)

include_directories(${EIGEN3_INCLUDE_DIRS})
message(STATUS "EIGEN3_INCLUDE_DIRS: ${EIGEN3_INCLUDE_DIRS}")

include_directories("include")
include_directories("include/Geometry")
include_directories("include/Render")
include_directories("include/GUI")
include_directories("include/Audio")
include_directories("include/Utils")
include_directories("include/Modal")
include_directories("include/Audio")

include_directories("external")
include_directories("external/json")
# include_directories("external/imgui")
include_directories("external/imgui/backends")
include_directories("external/libigl/include")
include_directories("external/spectra/include")
include_directories("external/argparse")

# generated from https://glad.dav1d.de/
include_directories("external/glad")
include_directories("external/miniaudio")

# include_directories("external/fTetWild/src")
target_include_directories(miniaudio_libvorbis PUBLIC /opt/homebrew/include)

# Add source files
file(GLOB_RECURSE SRCS "src/*.cpp")
file(GLOB_RECURSE GEOMETRY_SRCS "src/Geometry/*.cpp")
file(GLOB_RECURSE UTILS_SRCS "src/Utils/*.cpp")
file(GLOB_RECURSE MODAL_SRCS "src/Modal/*.cpp")
file(GLOB_RECURSE AUDIO_SRCS "src/Audio/*.cpp")
file(GLOB_RECURSE TRANSFER_SRCS "src/Transfer/*.cpp")


set(GLAD_SOURCES external/glad/glad.c)


# add_executable(test_cache test/geometry/test_cache.cpp ${UTILS_SRCS}
#                           ${GEOMETRY_SRCS}  ${GLAD_SOURCES})
# target_link_libraries(test_cache PUBLIC igl::glfw igl::opengl igl::imgui
#                                         igl_copyleft::tetgen)

# add_executable(test_picking test/modal/test_picking.cpp ${GLAD_SOURCES} ${MODAL_SRCS})
# target_link_libraries(test_picking PUBLIC igl::glfw igl::opengl igl::imgui
#                                           igl::embree)

add_executable(test_realtime_modal_sound test/modal/test_realtime_modal_sound.cpp ${GEOMETRY_SRCS} ${UTILS_SRCS} ${MODAL_SRCS} ${AUDIO_SRCS} ${TRANSFER_SRCS})
target_link_libraries(test_realtime_modal_sound PUBLIC igl::glfw igl::opengl igl::imgui igl::embree igl_copyleft::tetgen miniaudio)
# ---------------------------------------------------

target_include_directories(test_realtime_modal_sound PRIVATE /opt/homebrew/include)

# add_executable(test_miniaudio_white_noise test/miniaudio/procedural_whitenoise.cpp)
# target_link_libraries(test_miniaudio_white_noise PUBLIC miniaudio)

# add_executable(test_miniaudio_damped_sin_wave test/miniaudio/procedural_dampedSinWave.cpp)
# target_link_libraries(test_miniaudio_damped_sin_wave PUBLIC miniaudio)
