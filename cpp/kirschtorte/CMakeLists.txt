cmake_minimum_required(VERSION 3.18)

add_library(_kirschtorte_core MODULE
  kirschtorte.cpp
  decorators.cpp
  logger.cpp
  ilp_approximators.cpp
  decode_bitset.cpp
)

target_include_directories(_kirschtorte_core PRIVATE
  include
  "${boost_mirror_SOURCE_DIR}"
)

target_link_libraries(_kirschtorte_core PRIVATE
  Python3::Module
  ortools::ortools
  pybind11::module
)

if(NOT Python3_EXTENSION_SUFFIX)
    set(Python3_EXTENSION_SUFFIX ".so")
endif()

set_target_properties(_kirschtorte_core PROPERTIES
  PREFIX ""
  OUTPUT_NAME "_kirschtorte_core"
  SUFFIX "${Python3_EXTENSION_SUFFIX}"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)
