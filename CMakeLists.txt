# =============================================================================
# MyIndustrialCAD - Professional CAD System
# =============================================================================
cmake_minimum_required(VERSION 3.15)
project(MyIndustrialCAD
    VERSION 1.0.0
    DESCRIPTION "Professional CAD System with Boolean Operations"
    LANGUAGES CXX
)

# =============================================================================
# COMPILER SETTINGS
# =============================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Release mode optimizations
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(/O2 /Oi /GL)
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
    set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /LTCG")
    set(CMAKE_STATIC_LINKER_FLAGS_RELEASE "${CMAKE_STATIC_LINKER_FLAGS_RELEASE} /LTCG")
endif()

# Debug mode settings
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(/Zi /Od /RTC1)
    add_definitions(-D_DEBUG)
endif()

# Platform-specific settings
if(WIN32)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS NOMINMAX WIN32_LEAN_AND_MEAN)
    add_compile_options(/W4 /WX)  # High warning level, treat warnings as errors
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Werror)
endif()

# =============================================================================
# DEPENDENCIES
# =============================================================================
option(USE_CGAL "Use CGAL library for advanced Boolean operations" OFF)
option(BUILD_TESTS "Build test executables" ON)
option(BUILD_EXAMPLES "Build example programs" ON)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" OFF)
option(BUILD_GUI "Build GUI components" OFF)

# CGAL
if(USE_CGAL)
    find_package(CGAL REQUIRED COMPONENTS Core)
    if(CGAL_FOUND)
        add_definitions(-DHAVE_CGAL=1)
        message(STATUS "CGAL found: ${CGAL_VERSION}")
        include(${CGAL_USE_FILE})
    else()
        message(WARNING "CGAL not found - using internal Boolean operations")
        set(USE_CGAL OFF)
    endif()
else()
    message(STATUS "CGAL disabled - using internal Boolean operations")
endif()

# OpenMP (for parallel processing)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    add_definitions(-DUSE_OPENMP=1)
    message(STATUS "OpenMP found - enabling parallel processing")
else()
    message(STATUS "OpenMP not found - using single-threaded processing")
endif()

# =============================================================================
# PROJECT STRUCTURE
# =============================================================================
# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/core
    ${CMAKE_SOURCE_DIR}/core/part
    ${CMAKE_SOURCE_DIR}/core/sketcher
    ${CMAKE_SOURCE_DIR}/core/assembly
    ${CMAKE_SOURCE_DIR}/core/io
    ${CMAKE_SOURCE_DIR}/core/part_design
    ${CMAKE_SOURCE_DIR}/core/sheet_metal
    ${CMAKE_SOURCE_DIR}/core/surfacing
    ${CMAKE_SOURCE_DIR}/ui
    ${CMAKE_SOURCE_DIR}/python_api
)

# Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Create output directories
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
file(MAKE_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

# =============================================================================
# SUBDIRECTORIES
# =============================================================================
add_subdirectory(core)
add_subdirectory(tests)
if(EXISTS ${CMAKE_SOURCE_DIR}/ui/CMakeLists.txt)
    add_subdirectory(ui)
endif()

# =============================================================================
# MAIN EXECUTABLES
# =============================================================================
# Main CAD application (if we had one)
# add_executable(myindustrialcad main.cpp)
# target_link_libraries(myindustrialcad core_part core_sketcher)

# =============================================================================
# PACKAGE CONFIGURATION
# =============================================================================
# Installation
install(TARGETS core_part core_sketcher
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY core/ DESTINATION include/myindustrialcad
    FILES_MATCHING PATTERN "*.h"
)

# CPack configuration
set(CPACK_PACKAGE_NAME "MyIndustrialCAD")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_DESCRIPTION "Professional CAD System")
set(CPACK_PACKAGE_VENDOR "MyIndustrialCAD Team")
if(EXISTS ${CMAKE_SOURCE_DIR}/LICENSE.txt)
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE.txt")
else()
    # Temporäre Lizenzdatei erstellen
    file(WRITE ${CMAKE_BINARY_DIR}/LICENSE.txt "MyIndustrialCAD - Proprietary Software\nCopyright (c) 2025\n")
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_BINARY_DIR}/LICENSE.txt")
endif()
set(CPACK_NSIS_MODIFY_PATH ON)

if(WIN32)
    set(CPACK_GENERATOR "NSIS;ZIP")
else()
    set(CPACK_GENERATOR "TGZ")
endif()

include(CPack)

# =============================================================================
# CUSTOM TARGETS
# =============================================================================
# Documentation target
add_custom_target(docs
    COMMAND echo "Generating documentation..."
    COMMENT "Building documentation"
)

# Clean target
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_SOURCE_DIR}/out
    COMMENT "Cleaning all build artifacts"
)

# Format code target
find_program(CLANG_FORMAT_EXE NAMES clang-format clang-format-12 clang-format-11)
if(CLANG_FORMAT_EXE)
    add_custom_target(format
        COMMAND ${CLANG_FORMAT_EXE} -i --style=file ${CMAKE_SOURCE_DIR}/core/*.cpp ${CMAKE_SOURCE_DIR}/core/*/*.cpp ${CMAKE_SOURCE_DIR}/core/*/*.h
        COMMAND ${CLANG_FORMAT_EXE} -i --style=file ${CMAKE_SOURCE_DIR}/tests/*.cpp
        COMMENT "Formatting source code"
    )
endif()

# =============================================================================
# FINAL MESSAGE
# =============================================================================
message(STATUS "")
message(STATUS "===============================================")
message(STATUS "MyIndustrialCAD Configuration Summary")
message(STATUS "===============================================")
message(STATUS "Build type:          ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix:      ${CMAKE_INSTALL_PREFIX}")
message(STATUS "C++ standard:        ${CMAKE_CXX_STANDARD}")
message(STATUS "CGAL support:        ${USE_CGAL}")
message(STATUS "OpenMP support:      ${OpenMP_CXX_FOUND}")
message(STATUS "Build tests:         ${BUILD_TESTS}")
message(STATUS "Build examples:      ${BUILD_EXAMPLES}")
message(STATUS "===============================================")