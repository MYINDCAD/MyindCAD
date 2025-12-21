# C:\MyIndustrialCAD\CMakeLists.txt - KORRIGIERT
cmake_minimum_required(VERSION 3.15)
project(MyIndustrialCAD)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build-Optionen
option(BUILD_TESTS "Build tests" OFF)
option(BUILD_UI "Build UI components" OFF)
option(BUILD_SKETCHER "Build sketcher module" OFF)  # AUSGESCHALTET

# Warnung für CGAL
message(WARNING "CGAL not found - Boolean operations will use simplified algorithms")

# Setze Ausgabeverzeichnis
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/out)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/out)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/out)

# Unterverzeichnisse - NUR core/part
add_subdirectory(core/part)

# Sketcher AUSGESCHALTET
# if(BUILD_SKETCHER)
#     add_subdirectory(sketcher)
# endif()

if(BUILD_UI)
    add_subdirectory(ui)
endif()

if(BUILD_TESTS)
    add_subdirectory(tests)
endif()


# Tests
add_subdirectory(tests)  # DIESE ZEILE MUSS VORHANDEN SEIN!