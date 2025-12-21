#include "../core/part/Part.h"
#include <iostream>

void testBooleanUnion() {
    std::cout << "=== Test Boolean Union ===\n";

    // Zwei überlappende Boxen erstellen
    auto box1 = Part::createBox(10, 10, 10);
    auto box2 = Part::createBox(5, 5, 15);  // Überschneidung

    std::cout << "Box1: 10x10x10\n";
    std::cout << "Box2: 5x5x15\n";
    std::cout << "Performing union...\n";

    try {
        auto unionSolid = Part::booleanUnion(*box1, *box2);

        std::cout << "Union successful!\n";
        std::cout << "Result solid: " << unionSolid->getName() << "\n";

        Point3D min, max;
        unionSolid->getBoundingBox(min, max);
        std::cout << "Bounding Box: [" << min.x << "," << min.y << "," << min.z
            << "] to [" << max.x << "," << max.y << "," << max.z << "]\n";

    }
    catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
    }
}

void testBooleanDifference() {
    std::cout << "\n=== Test Boolean Difference ===\n";

    auto box1 = Part::createBox(10, 10, 10);
    auto box2 = Part::createBox(5, 5, 5);

    std::cout << "Box1 - Box2 (subtract smaller from larger)\n";

    try {
        auto diffSolid = Part::booleanDifference(*box1, *box2);

        std::cout << "Difference successful!\n";
        std::cout << "Result solid: " << diffSolid->getName() << "\n";

    }
    catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
    }
}

void testBooleanIntersection() {
    std::cout << "\n=== Test Boolean Intersection ===\n";

    auto box1 = Part::createBox(10, 10, 10);
    auto box2 = Part::createBox(5, 5, 15);  // Teilweise Überlappung

    std::cout << "Box1 ∩ Box2 (intersection)\n";

    try {
        auto interSolid = Part::booleanIntersection(*box1, *box2);

        std::cout << "Intersection successful!\n";
        std::cout << "Result solid: " << interSolid->getName() << "\n";

        if (!interSolid->getVertices().empty()) {
            std::cout << "Intersection volume: " << interSolid->volume() << " mm³\n";
        }

    }
    catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
    }
}

int main() {
    std::cout << "=======================================\n";
    std::cout << "MyIndustrialCAD - Boolean Operations Tests\n";
    std::cout << "=======================================\n\n";

    testBooleanUnion();
    testBooleanDifference();
    testBooleanIntersection();

    std::cout << "\n=== All boolean tests completed ===\n";
    return 0;
}