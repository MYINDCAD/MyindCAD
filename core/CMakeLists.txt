# =============================================================================
# Core Module - MyIndustrialCAD
# =============================================================================
cmake_minimum_required(VERSION 3.15)

# Core library version
set(CORE_VERSION_MAJOR 1)
set(CORE_VERSION_MINOR 0)
set(CORE_VERSION_PATCH 0)
set(CORE_VERSION ${CORE_VERSION_MAJOR}.${CORE_VERSION_MINOR}.${CORE_VERSION_PATCH})

# Export symbols for shared library
if(BUILD_SHARED_LIBS)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
    set(CMAKE_CXX_VISIBILITY_PRESET hidden)
    set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)
endif()

# =============================================================================
# SUBDIRECTORIES
# =============================================================================
# Add all core modules
add_subdirectory(part)

# Optional modules
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/sketcher/CMakeLists.txt)
    add_subdirectory(sketcher)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/assembly/CMakeLists.txt)
    add_subdirectory(assembly)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/io/CMakeLists.txt)
    add_subdirectory(io)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/part_design/CMakeLists.txt)
    add_subdirectory(part_design)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/sheet_metal/CMakeLists.txt)
    add_subdirectory(sheet_metal)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/surfacing/CMakeLists.txt)
    add_subdirectory(surfacing)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Python_api/CMakeLists.txt)
    add_subdirectory(Python_api)
endif()

# =============================================================================
# CORE CONFIGURATION
# =============================================================================
# Common compile definitions for all core modules
add_compile_definitions(
    MYINDUSTRIALCAD_CORE
    CORE_VERSION=${CORE_VERSION}
)

# Common include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/part
    ${CMAKE_CURRENT_SOURCE_DIR}/sketcher
)

# =============================================================================
# EXPORT TARGETS
# =============================================================================
# Create an interface library for easy linking
add_library(core INTERFACE)
add_library(MyIndustrialCAD::core ALIAS core)

target_include_directories(core INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

# Dependencies
if(USE_CGAL)
    target_link_libraries(core INTERFACE CGAL::CGAL)
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(core INTERFACE OpenMP::OpenMP_CXX)
endif()

# =============================================================================
# INSTALLATION
# =============================================================================
# Export targets for use from build tree
export(TARGETS core
    FILE "${CMAKE_BINARY_DIR}/MyIndustrialCADCoreTargets.cmake"
)

# Install configuration
install(TARGETS core
    EXPORT MyIndustrialCADCoreTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(EXPORT MyIndustrialCADCoreTargets
    FILE MyIndustrialCADCoreConfig.cmake
    DESTINATION lib/cmake/MyIndustrialCAD
)

# Install headers
install(DIRECTORY part/
    DESTINATION include/myindustrialcad/part
    FILES_MATCHING PATTERN "*.h"
)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/sketcher)
    install(DIRECTORY sketcher/
        DESTINATION include/myindustrialcad/sketcher
        FILES_MATCHING PATTERN "*.h"
    )
endif()