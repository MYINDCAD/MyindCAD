# Library for the Sketcher module
add_library(core_sketcher
    Sketch.cpp
)

# Include path for headers
target_include_directories(core_sketcher PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
