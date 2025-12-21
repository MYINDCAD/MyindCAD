# C:\MyIndustrialCAD\core\part\CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(core_part)

# Header-Dateien
set(HEADER_FILES
    Part.h
    BooleanOperations.h
    Solid.h
    Vertex.h
    Edge.h
    Face.h
)

# Quell-Dateien - Solid.cpp entfernen da Solid.h inline ist
set(SOURCE_FILES
    Part.cpp
    BooleanOperations.cpp
    # Solid.cpp NICHT MEHR - alles in Solid.h inline
)

# Bibliothek erstellen
add_library(core_part STATIC ${SOURCE_FILES})

# Include-Verzeichnisse
target_include_directories(core_part
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/../sketcher
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/../..
)

# Compiler-Flags
target_compile_definitions(core_part PRIVATE _USE_MATH_DEFINES)

# C++ Standard setzen
set_target_properties(core_part PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

# Bibliothek exportieren (optional)
set(CORE_PART_LIB core_part PARENT_SCOPE)