# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/daiki/Documents/Stanford/Research/SimpleModal/build/_deps/imguizmo-src"
  "/Users/daiki/Documents/Stanford/Research/SimpleModal/build/_deps/imguizmo-build"
  "/Users/daiki/Documents/Stanford/Research/SimpleModal/build/_deps/imguizmo-subbuild/imguizmo-populate-prefix"
  "/Users/daiki/Documents/Stanford/Research/SimpleModal/build/_deps/imguizmo-subbuild/imguizmo-populate-prefix/tmp"
  "/Users/daiki/Documents/Stanford/Research/SimpleModal/build/_deps/imguizmo-subbuild/imguizmo-populate-prefix/src/imguizmo-populate-stamp"
  "/Users/daiki/Documents/Stanford/Research/SimpleModal/build/_deps/imguizmo-subbuild/imguizmo-populate-prefix/src"
  "/Users/daiki/Documents/Stanford/Research/SimpleModal/build/_deps/imguizmo-subbuild/imguizmo-populate-prefix/src/imguizmo-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/daiki/Documents/Stanford/Research/SimpleModal/build/_deps/imguizmo-subbuild/imguizmo-populate-prefix/src/imguizmo-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/daiki/Documents/Stanford/Research/SimpleModal/build/_deps/imguizmo-subbuild/imguizmo-populate-prefix/src/imguizmo-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
