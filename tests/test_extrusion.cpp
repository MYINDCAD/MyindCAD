#include "../core/part/Part.h"
#include "../core/sketcher/Sketch.h"
#include <iostream>

void testRectangleExtrusion() {
    std::cout << "=== Test 1: Rectangle Extrusion ===\n";

    Sketch sketch;
    sketch.addRectangle(0, 0, 10, 10);

    std::cout << "Sketch: 10x10 rectangle\n";
    std::cout << "Area: " << sketch.calculateArea() << " mm²\n";
    std::cout << "Extruding to height 5mm...\n";

    try {
        auto solid = Part::extrudeSketch(sketch, 5.0);

        std::cout << "✓ Extrusion successful!\n";
        std::cout << "Solid name: " << solid->getName() << "\n";
        std::cout << "Faces: " << solid->getFaces().size() << " (expected: 6)\n";
        std::cout << "Vertices: " << solid->getVertices().size() << " (expected: 8)\n";
        std::cout << "Edges: " << solid->getEdges().size() << " (expected: 12)\n";
        std::cout << "Volume: " << solid->volume() << " mm³ (expected: 500)\n";

        // Bounding Box prüfen
        Point3D min, max;
        solid->getBoundingBox(min, max);
        std::cout << "Bounding Box: [" << min.x << "," << min.y << "," << min.z
            << "] to [" << max.x << "," << max.y << "," << max.z << "]\n";

    }
    catch (const std::exception& e) {
        std::cout << "✗ Error: " << e.what() << "\n";
    }
}

void testCircleExtrusion() {
    std::cout << "\n=== Test 2: Circle Extrusion ===\n";

    Sketch sketch;
    sketch.addCircle(0, 0, 5.0);  // Radius 5

    std::cout << "Sketch: Circle radius 5\n";
    std::cout << "Area: " << sketch.calculateArea() << " mm² (expected: ~78.54)\n";
    std::cout << "Extruding to height 10mm...\n";

    try {
        auto solid = Part::extrudeSketch(sketch, 10.0);

        std::cout << "✓ Extrusion successful!\n";
        std::cout << "Solid name: " << solid->getName() << "\n";
        // Kreisextrusion hat viele Vertices (aufgrund der Diskretisierung)
        std::cout << "Faces: " << solid->getFaces().size() << "\n";
        std::cout << "Volume: " << solid->volume() << " mm³\n";

    }
    catch (const std::exception& e) {
        std::cout << "✗ Error: " << e.what() << "\n";
    }
}

void testInvalidExtrusion() {
    std::cout << "\n=== Test 3: Invalid Sketch ===\n";

    Sketch sketch;
    sketch.addLine(0, 0, 10, 0);  // Nur eine Linie - nicht geschlossen

    std::cout << "Sketch: Single line (not closed)\n";
    std::cout << "Is closed: " << (sketch.isClosed() ? "Yes" : "No") << "\n";
    std::cout << "Attempting extrusion...\n";

    try {
        auto solid = Part::extrudeSketch(sketch, 5.0);
        std::cout << "✗ ERROR: Should have thrown exception!\n";
    }
    catch (const std::exception& e) {
        std::cout << "✓ Correctly rejected: " << e.what() << "\n";
    }
}

int main() {
    std::cout << "========================================\n";
    std::cout << "MyIndustrialCAD - Extrusion Module Tests\n";
    std::cout << "========================================\n\n";

    testRectangleExtrusion();
    testCircleExtrusion();
    testInvalidExtrusion();

    std::cout << "\n=== All extrusion tests completed ===\n";
    return 0;
}