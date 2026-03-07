cmake_minimum_required(VERSION 3.18)
project(allegro LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-Os -ffunction-sections -fdata-sections)
endif()

find_package(Threads REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(ZLIB REQUIRED)

include(FetchContent)

# =========================
# pybind11 (REQUIRED for kirschtorte glue code)
# =========================
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v3.0.2
)
set(PYBIND11_FINDPYTHON ON)
FetchContent_MakeAvailable(pybind11)

# =========================
# OR-Tools
# =========================
set(BUILD_DEPS ON CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_SAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

# Solver settings
set(USE_COINOR OFF CACHE BOOL "" FORCE)
set(USE_HIGHS OFF CACHE BOOL "" FORCE)
set(USE_SCIP OFF CACHE BOOL "" FORCE)
set(USE_BOP OFF CACHE BOOL "" FORCE)
set(USE_MATH_OPT OFF CACHE BOOL "" FORCE)
set(USE_PDLP OFF CACHE BOOL "" FORCE)
set(USE_DOTNET_8 OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
  or-tools
  URL https://github.com/google/or-tools/archive/refs/tags/v9.11.tar.gz
)
FetchContent_MakeAvailable(or-tools)

# -------------------------
# 1) kirschtorte Python extension
# -------------------------
add_subdirectory(kirschtorte) 

# -------------------------
# 2) RNAstructure
# -------------------------
include(ExternalProject)
set(RNASTRUCTURE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/uCRISPR/RNAstructure")
set(RNASTRUCTURE_LIB "${CMAKE_BINARY_DIR}/librnastructure.a")

ExternalProject_Add(
  rnastructure_build
  SOURCE_DIR "${RNASTRUCTURE_DIR}"
  CONFIGURE_COMMAND ""
  BUILD_IN_SOURCE 1
  BUILD_COMMAND ${CMAKE_COMMAND} -E env CFLAGS=-fPIC CXXFLAGS=-fPIC make scorer
  INSTALL_COMMAND ""
)

ExternalProject_Add_Step(
  rnastructure_build archive
  DEPENDEES build
  COMMAND ${CMAKE_COMMAND} -E rm -f "${RNASTRUCTURE_LIB}"
  COMMAND ${CMAKE_COMMAND} -E env bash -c
    "cd \"${RNASTRUCTURE_DIR}\" && find . -type f -name '*.o' -not -path './exe/*' | xargs ar rcs \"${RNASTRUCTURE_LIB}\" && ranlib \"${RNASTRUCTURE_LIB}\""
)

add_library(rnastructure STATIC IMPORTED GLOBAL)
set_target_properties(rnastructure PROPERTIES IMPORTED_LOCATION "${RNASTRUCTURE_LIB}")
add_dependencies(rnastructure rnastructure_build)

# -------------------------
# 3) uCRISPR_scorer
# -------------------------
add_executable(uCRISPR_scorer
  "uCRISPR/main.cpp"
  "uCRISPR/uCRISPR_scorer.cpp"
)

target_include_directories(uCRISPR_scorer PRIVATE
  "uCRISPR"
  "uCRISPR/include"
  "uCRISPR/RNAstructure"
  "boost_1_81_0"
)

target_link_libraries(uCRISPR_scorer PRIVATE rnastructure ZLIB::ZLIB Threads::Threads)

set_target_properties(uCRISPR_scorer PROPERTIES 
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)
