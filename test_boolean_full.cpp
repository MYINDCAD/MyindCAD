// C:\MyIndustrialCAD\test_boolean_full.cpp
#include <iostream>
#include <memory>

// Includes für unsere Bibliotheken
#include "core/part/Part.h"
#include "core/part/Solid.h"
#include "core/part/BooleanOperations.h"

int main() {
    std::cout << "===============================================\n";
    std::cout << "    BOOLEAN OPERATIONS TEST - PHASE 2\n";
    std::cout << "===============================================\n\n";

    try {
        std::cout << "1. Creating test solids...\n";

        // Erstelle zwei Test-Solids (Boxen)
        auto box1 = Part::createBox(10, 10, 10);
        auto box2 = Part::createBox(10, 10, 10);

        // Box2 verschieben, damit sie sich überlappen
        box2->translate(Vector3D(5, 5, 5));

        std::cout << "   Box1 created: " << box1->getName() << "\n";
        std::cout << "   Box2 created: " << box2->getName() << "\n";
        std::cout << "   Box2 translated by (5, 5, 5)\n\n";

        std::cout << "2. Testing Boolean Operations...\n";

        // Test 1: Union
        std::cout << "\n   a) Testing UNION operation:\n";
        std::cout << "   ---------------------------\n";
        auto unionResult = BooleanOperations::booleanUnion(*box1, *box2);
        if (unionResult) {
            std::cout << "   ✓ Union created: " << unionResult->getName() << "\n";
            std::cout << "   Vertices: " << unionResult->getVertices().size() << "\n";
            std::cout << "   Edges: " << unionResult->getEdges().size() << "\n";
            std::cout << "   Faces: " << unionResult->getFaces().size() << "\n";
        }
        else {
            std::cout << "   ✗ Union failed!\n";
        }

        // Test 2: Difference
        std::cout << "\n   b) Testing DIFFERENCE operation:\n";
        std::cout << "   --------------------------------\n";
        auto diffResult = BooleanOperations::booleanDifference(*box1, *box2);
        if (diffResult) {
            std::cout << "   ✓ Difference created: " << diffResult->getName() << "\n";
            std::cout << "   Vertices: " << diffResult->getVertices().size() << "\n";
            std::cout << "   Edges: " << diffResult->getEdges().size() << "\n";
            std::cout << "   Faces: " << diffResult->getFaces().size() << "\n";
        }
        else {
            std::cout << "   ✗ Difference failed!\n";
        }

        // Test 3: Intersection
        std::cout << "\n   c) Testing INTERSECTION operation:\n";
        std::cout << "   -----------------------------------\n";
        auto intersectResult = BooleanOperations::booleanIntersection(*box1, *box2);
        if (intersectResult) {
            std::cout << "   ✓ Intersection created: " << intersectResult->getName() << "\n";
            std::cout << "   Vertices: " << intersectResult->getVertices().size() << "\n";
            std::cout << "   Edges: " << intersectResult->getEdges().size() << "\n";
            std::cout << "   Faces: " << intersectResult->getFaces().size() << "\n";
        }
        else {
            std::cout << "   ✗ Intersection failed!\n";
        }

        std::cout << "\n3. Testing via Part class wrapper...\n";

        // Test via Part-Klasse (Wrapper-Funktionen)
        Part part("TestPart");
        part.addSolid(box1);

        auto unionViaPart = Part::booleanUnion(*box1, *box2);
        if (unionViaPart) {
            std::cout << "   ✓ Part::booleanUnion() works\n";
        }

        std::cout << "\n4. Solid properties test...\n";

        // Teste Solid-Eigenschaften
        Point3D min, max;
        box1->getBoundingBox(min, max);
        std::cout << "   Box1 Bounding Box:\n";
        std::cout << "     Min: (" << min.x << ", " << min.y << ", " << min.z << ")\n";
        std::cout << "     Max: (" << max.x << ", " << max.y << ", " << max.z << ")\n";

        std::cout << "\n   Box1 Volume: " << box1->volume() << " mm³\n";
        std::cout << "   Box1 Valid: " << (box1->isValid() ? "Yes" : "No") << "\n";

        std::cout << "\n===============================================\n";
        std::cout << "    ✅ BOOLEAN OPERATIONS TEST COMPLETED\n";
        std::cout << "===============================================\n";

    }
    catch (const std::exception& e) {
        std::cout << "\n❌ ERROR: " << e.what() << "\n";
        return 1;
    }
    catch (...) {
        std::cout << "\n❌ Unknown error occurred\n";
        return 1;
    }

    return 0;
}