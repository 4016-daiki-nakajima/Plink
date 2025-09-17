cmake_minimum_required(VERSION 3.10)
project(Plink)

# -----------------------------------------------------------------------------
# General 
# -----------------------------------------------------------------------------

set(PROJECT_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
add_definitions(-DPROJECT_ROOT_DIR="${PROJECT_ROOT_DIR}")

set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS TRUE)
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------

# Require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Core dependencies
find_package(OpenGL REQUIRED)
find_package(Eigen3 REQUIRED)
# find_package(Boost REQUIRED)

# Header-only libraries
add_subdirectory(external/libigl)
igl_include(glfw)
igl_include(opengl)
igl_include(imgui)
igl_include(embree)
igl_include(copyleft tetgen)

add_subdirectory(external/spectra)
add_subdirectory(external/miniaudio)

include_directories("include")

include_directories("external")
include_directories("external/json")
include_directories("external/spectra/include")
include_directories("external/argparse")
include_directories("external/glad")
include_directories("external/miniaudio")

# -----------------------------------------------------------------------------
# Build main executable
# -----------------------------------------------------------------------------

file(GLOB_RECURSE GEOMETRY_SRCS "src/Geometry/*.cpp")
file(GLOB_RECURSE UTILS_SRCS "src/Utils/*.cpp")
file(GLOB_RECURSE MODAL_SRCS "src/Modal/*.cpp")
file(GLOB_RECURSE AUDIO_SRCS "src/Audio/*.cpp")
file(GLOB_RECURSE TRANSFER_SRCS "src/Transfer/*.cpp")
set(GLAD_SOURCES external/glad/glad.c)

add_executable(Plink 
    main.cpp
    ${GEOMETRY_SRCS}
    ${UTILS_SRCS}
    ${MODAL_SRCS}
    ${AUDIO_SRCS}
    ${TRANSFER_SRCS}
    ${GLAD_SOURCES}
)

target_link_libraries(Plink PUBLIC
    igl::glfw
    igl::opengl
    igl::imgui
    igl::embree
    igl_copyleft::tetgen
    miniaudio
)
