# Install script for directory: /Users/daiki/Documents/Stanford/Research/SimpleModal/build/_deps/embree-src/common/lexers

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Library/Developer/CommandLineTools/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/daiki/Documents/Stanford/Research/SimpleModal/build/_deps/embree-build/liblexers.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/liblexers.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/liblexers.a")
    execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/liblexers.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "devel" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/embree-4.4.0/lexers-targets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/embree-4.4.0/lexers-targets.cmake"
         "/Users/daiki/Documents/Stanford/Research/SimpleModal/build/_deps/embree-build/common/lexers/CMakeFiles/Export/bcb401d0aa29ca560a1d675a6e229d39/lexers-targets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/embree-4.4.0/lexers-targets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/embree-4.4.0/lexers-targets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/embree-4.4.0" TYPE FILE FILES "/Users/daiki/Documents/Stanford/Research/SimpleModal/build/_deps/embree-build/common/lexers/CMakeFiles/Export/bcb401d0aa29ca560a1d675a6e229d39/lexers-targets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/embree-4.4.0" TYPE FILE FILES "/Users/daiki/Documents/Stanford/Research/SimpleModal/build/_deps/embree-build/common/lexers/CMakeFiles/Export/bcb401d0aa29ca560a1d675a6e229d39/lexers-targets-release.cmake")
  endif()
endif()

