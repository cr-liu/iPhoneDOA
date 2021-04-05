cmake_minimum_required(VERSION 3.5)

project(iphone-doa LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

include_directories(/usr/local/include)
link_directories(/usr/local/lib)

add_library(loclib)
target_sources(loclib
  PRIVATE
  src/loclib/arraymanifold.hpp
  src/loclib/arraymanifold.cpp
  src/loclib/dsloc.hpp
  src/loclib/dsloc.cpp
  src/loclib/musicloc.hpp
  src/loclib/musicloc.cpp
)
target_link_libraries(loclib sndfile)

add_executable(iphone-doa
  src/main.cpp
  src/wavfile.hpp
  src/wavfile.cpp
  src/winfunc.hpp
  src/winfunc.cpp
  src/corrmat.hpp
  src/corrmat.cpp
)

target_link_libraries(iphone-doa sndfile)
target_link_libraries(iphone-doa fftw3f)
target_link_libraries(iphone-doa loclib)
