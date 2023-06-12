cmake_minimum_required(VERSION 3.18)

option(CMAKE_EXPORT_COMPILE_COMMANDS "Export compile command" OFF)

project(kirschtorte VERSION 1.0 LANGUAGES CXX)
message(STATUS "${PROJECT_NAME} version: ${PROJECT_VERSION}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True)

# Load OR-Tools
find_package(ortools REQUIRED CONFIG)

# Load a Python interpreter
find_package(PythonInterp REQUIRED)

# Directory cointaining UseCython.cmake and FindCython.cmake
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/cmake/cmake-cython)

# This defines cython_add_module
include(UseCython)

# Tells UseCython to compile this file as a c++ file
set_source_files_properties(kirschtorte.pyx PROPERTIES CYTHON_IS_CXX TRUE)

# Create python module
cython_add_module(${PROJECT_NAME} kirschtorte.pyx kirschtorte.cpp decorators.cpp decode_bitset.cpp logger.cpp ilp_approximators.cpp)

# Location of .h files
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Link OR-Tools libs
target_link_libraries(${PROJECT_NAME} ortools::ortools)

# Boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${Boost_DYNAMIC_BITSET_LIBRARY})