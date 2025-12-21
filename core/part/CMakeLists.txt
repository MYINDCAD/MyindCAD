# =============================================================================
# Part Module - Solid Modeling and Boolean Operations
# =============================================================================
cmake_minimum_required(VERSION 3.15)

# Module info
set(PART_MODULE_NAME core_part)
set(PART_MODULE_DESCRIPTION "Solid Modeling and Boolean Operations Module")
set(PART_MODULE_VERSION 1.0.0)

# =============================================================================
# SOURCE FILES
# =============================================================================
set(PART_SOURCES
    Part.cpp
    BooleanOperations.cpp
)

set(PART_HEADERS
    Part.h
    BooleanOperations.h
    Solid.h
    Face.h
    Edge.h
    Vertex.h
    Point3D.h
)

# Optional: STL Export if file exists
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/STLExport.cpp)
    list(APPEND PART_SOURCES STLExport.cpp)
    list(APPEND PART_HEADERS STLExport.h)
    add_definitions(-DHAVE_STL_EXPORT=1)
endif()

# Optional: SketcherToSolid if file exists
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/SketcherToSolid.cpp)
    list(APPEND PART_SOURCES SketcherToSolid.cpp)
    list(APPEND PART_HEADERS SketcherToSolid.h)
    add_definitions(-DHAVE_SKETCHER_TO_SOLID=1)
endif()

# =============================================================================
# LIBRARY TARGET
# =============================================================================
# Create library
add_library(${PART_MODULE_NAME} STATIC ${PART_SOURCES})

# Set library properties
set_target_properties(${PART_MODULE_NAME} PROPERTIES
    VERSION ${PART_MODULE_VERSION}
    SOVERSION 1
    OUTPUT_NAME "core_part"
    ARCHIVE_OUTPUT_NAME "core_part"
    LIBRARY_OUTPUT_NAME "core_part"
    RUNTIME_OUTPUT_NAME "core_part"
    DEBUG_POSTFIX "d"
    RELEASE_POSTFIX ""
    MINSIZEREL_POSTFIX "minsize"
    RELWITHDEBINFO_POSTFIX "relwithdebinfo"
)

# Alias for consistent naming
add_library(MyIndustrialCAD::core_part ALIAS ${PART_MODULE_NAME})

# =============================================================================
# INCLUDE DIRECTORIES
# =============================================================================
target_include_directories(${PART_MODULE_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/myindustrialcad/part>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# =============================================================================
# COMPILE DEFINITIONS
# =============================================================================
target_compile_definitions(${PART_MODULE_NAME}
    PUBLIC
        PART_MODULE
        PART_VERSION=${PART_MODULE_VERSION}
    PRIVATE
        _USE_MATH_DEFINES
)

# =============================================================================
# COMPILE OPTIONS
# =============================================================================
if(MSVC)
    target_compile_options(${PART_MODULE_NAME}
        PRIVATE
            /W4     # High warning level
            /wd4251 # Disable DLL interface warning
            /wd4275 # Disable DLL interface warning
    )
else()
    target_compile_options(${PART_MODULE_NAME}
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wno-unused-parameter
    )
endif()

# =============================================================================
# LINK LIBRARIES
# =============================================================================
# Link with math library
if(UNIX AND NOT APPLE)
    target_link_libraries(${PART_MODULE_NAME} PRIVATE m)
endif()

# Link with CGAL if available
if(USE_CGAL)
    target_link_libraries(${PART_MODULE_NAME} PUBLIC CGAL::CGAL)
    target_compile_definitions(${PART_MODULE_NAME} PUBLIC USE_CGAL=1)
endif()

# Link with OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PART_MODULE_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

# =============================================================================
# INSTALLATION
# =============================================================================
# Install library
install(TARGETS ${PART_MODULE_NAME}
    EXPORT MyIndustrialCADPartTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

# Install headers
foreach(header ${PART_HEADERS})
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${header})
        install(FILES ${header}
            DESTINATION include/myindustrialcad/part
        )
    endif()
endforeach()

# Export targets
install(EXPORT MyIndustrialCADPartTargets
    FILE MyIndustrialCADPartConfig.cmake
    NAMESPACE MyIndustrialCAD::
    DESTINATION lib/cmake/MyIndustrialCAD
)

# =============================================================================
# TESTING
# =============================================================================
if(BUILD_TESTS)
    # Unit tests for Part module
    add_executable(test_part_module
        test.cpp  # Rename your test.cpp to something specific
    )
    
    target_link_libraries(test_part_module ${PART_MODULE_NAME})
    
    # Set output directory
    set_target_properties(test_part_module PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests
    )
    
    # Add test to CTest
    add_test(NAME test_part_module
        COMMAND test_part_module
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
    )
endif()

# =============================================================================
# EXPORT FOR PARENT
# =============================================================================
# Make module available to parent scope
set(PART_MODULE_TARGET ${PART_MODULE_NAME} PARENT_SCOPE)
set(PART_MODULE_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

# =============================================================================
# MODULE DOCUMENTATION
# =============================================================================
message(STATUS "Configured Part Module: ${PART_MODULE_NAME}")
message(STATUS "  Version: ${PART_MODULE_VERSION}")
message(STATUS "  Sources: ${PART_SOURCES}")
message(STATUS "  Headers: ${PART_HEADERS}")