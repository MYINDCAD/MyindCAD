# =============================================================================
# Tests for MyIndustrialCAD
# =============================================================================
cmake_minimum_required(VERSION 3.15)

# Only build tests if requested
if(NOT BUILD_TESTS)
    message(STATUS "Tests disabled (BUILD_TESTS=OFF)")
    return()
endif()

# Enable testing
enable_testing()

# Include CTest module
include(CTest)

# =============================================================================
# TEST SOURCE FILES
# =============================================================================
set(TEST_SOURCES
    test_all_operations.cpp
    test_boolean_full.cpp
)

# =============================================================================
# TEST CONFIGURATION
# =============================================================================
# Output directory for tests
set(TEST_OUTPUT_DIR ${CMAKE_BINARY_DIR}/tests)
file(MAKE_DIRECTORY ${TEST_OUTPUT_DIR})

# =============================================================================
# INDIVIDUAL TEST EXECUTABLES
# =============================================================================

# Test 1: Complete Boolean Operations Test
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test_all_operations.cpp)
    add_executable(test_all_operations
        test_all_operations.cpp
    )
    
    target_link_libraries(test_all_operations
        core_part
    )
    
    set_target_properties(test_all_operations PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${TEST_OUTPUT_DIR}
        OUTPUT_NAME test_all_operations
    )
    
    # Add to CTest
    add_test(NAME test_all_operations
        COMMAND test_all_operations
        WORKING_DIRECTORY ${TEST_OUTPUT_DIR}
    )
    
    # Set test properties
    set_tests_properties(test_all_operations PROPERTIES
        LABELS "boolean;operations;regression"
        TIMEOUT 30
    )
endif()

# Test 2: Boolean Full Test
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test_boolean_full.cpp)
    add_executable(test_boolean_full
        test_boolean_full.cpp
    )
    
    target_link_libraries(test_boolean_full
        core_part
    )
    
    set_target_properties(test_boolean_full PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${TEST_OUTPUT_DIR}
        OUTPUT_NAME test_boolean_full
    )
    
    add_test(NAME test_boolean_full
        COMMAND test_boolean_full
        WORKING_DIRECTORY ${TEST_OUTPUT_DIR}
    )
    
    set_tests_properties(test_boolean_full PROPERTIES
        LABELS "boolean;basic"
        TIMEOUT 15
    )
endif()

# Test 3: Phase 1 Success Test
if(EXISTS ${CMAKE_SOURCE_DIR}/phase1_success.cpp)
    add_executable(test_phase1_success
        ${CMAKE_SOURCE_DIR}/phase1_success.cpp
    )
    
    target_link_libraries(test_phase1_success
        core_part
    )
    
    set_target_properties(test_phase1_success PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${TEST_OUTPUT_DIR}
        OUTPUT_NAME test_phase1_success
    )
    
    add_test(NAME test_phase1_success
        COMMAND test_phase1_success
        WORKING_DIRECTORY ${TEST_OUTPUT_DIR}
    )
    
    set_tests_properties(test_phase1_success PROPERTIES
        LABELS "phase1;basic"
        TIMEOUT 10
    )
endif()

# =============================================================================
# TEST DEPENDENCIES
# =============================================================================
# Make sure core modules are built before tests
foreach(test_target test_all_operations test_boolean_full test_phase1_success)
    if(TARGET ${test_target})
        add_dependencies(${test_target} core_part)
    endif()
endforeach()

# =============================================================================
# TEST DISCOVERY (for future Google Test integration)
# =============================================================================
option(USE_GOOGLE_TEST "Use Google Test framework" OFF)

if(USE_GOOGLE_TEST)
    find_package(GTest REQUIRED)
    include(GoogleTest)
    
    # Example for future Google Test integration
    # add_executable(part_unit_tests part_unit_tests.cpp)
    # target_link_libraries(part_unit_tests core_part GTest::gtest GTest::gtest_main)
    # gtest_discover_tests(part_unit_tests)
endif()

# =============================================================================
# TEST RUNNER
# =============================================================================
# Custom target to run all tests
add_custom_target(run_tests
    COMMAND ctest --output-on-failure
    DEPENDS test_all_operations test_boolean_full test_phase1_success
    COMMENT "Running all tests"
)

# Custom target to run tests with verbose output
add_custom_target(run_tests_verbose
    COMMAND ctest --verbose
    DEPENDS test_all_operations test_boolean_full test_phase1_success
    COMMENT "Running all tests with verbose output"
)

# =============================================================================
# TEST DATA
# =============================================================================
# Copy test data to output directory
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test_data)
    add_custom_command(
        TARGET test_all_operations
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_CURRENT_SOURCE_DIR}/test_data
            ${TEST_OUTPUT_DIR}/test_data
        COMMENT "Copying test data"
    )
endif()

# =============================================================================
# TEST SUMMARY
# =============================================================================
message(STATUS "Configured tests:")
foreach(test_target test_all_operations test_boolean_full test_phase1_success)
    if(TARGET ${test_target})
        message(STATUS "  - ${test_target}")
    endif()
endforeach()