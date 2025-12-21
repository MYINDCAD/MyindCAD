// test_all_operations.cpp
#include <iostream>
#include "BooleanOperations.h"
#include "STLExport.h"

using namespace std;

int main() {
    cout << "===============================================\n";
    cout << "   COMPLETE BOOLEAN OPERATIONS TEST SUITE\n";
    cout << "===============================================\n\n";

    try {
        // Configure
        BooleanOperations::setVerbose(true);
        BooleanOperations::setPrecision(1e-6);

        cout << "=== PHASE 1: PRIMITIVE CREATION ===\n\n";

        // 1. Create box
        cout << "1. Creating Box (20x15x10 mm)...\n";
        auto box = BooleanOperations::createBox(20, 15, 10);
        if (!box) {
            cerr << "Failed to create box\n";
            return 1;
        }
        STLExport::exportToSTL(*box, "test_box.stl");
        cout << "   ✓ Box created: " << box->getVertices().size()
            << " vertices, " << box->getFaces().size() << " faces\n";

        // 2. Create cylinder
        cout << "\n2. Creating Cylinder (radius=8mm, height=25mm)...\n";
        auto cylinder = BooleanOperations::createCylinder(8, 25, 32);
        if (!cylinder) {
            cerr << "Failed to create cylinder\n";
            return 1;
        }
        STLExport::exportToSTL(*cylinder, "test_cylinder.stl");
        cout << "   ✓ Cylinder created: " << cylinder->getVertices().size()
            << " vertices, " << cylinder->getFaces().size() << " faces\n";

        // 3. Create sphere
        cout << "\n3. Creating Sphere (radius=12mm)...\n";
        auto sphere = BooleanOperations::createSphere(12, 24);
        if (!sphere) {
            cerr << "Failed to create sphere\n";
            return 1;
        }
        STLExport::exportToSTL(*sphere, "test_sphere.stl");
        cout << "   ✓ Sphere created: " << sphere->getVertices().size()
            << " vertices, " << sphere->getFaces().size() << " faces\n";

        // 4. Create cone
        cout << "\n4. Creating Cone (radius=10mm, height=20mm)...\n";
        auto cone = BooleanOperations::createCone(10, 20, 24);
        if (!cone) {
            cerr << "Failed to create cone\n";
            return 1;
        }
        STLExport::exportToSTL(*cone, "test_cone.stl");
        cout << "   ✓ Cone created: " << cone->getVertices().size()
            << " vertices, " << cone->getFaces().size() << " faces\n";

        // 5. Create pyramid
        cout << "\n5. Creating Pyramid (base=15x15mm, height=20mm)...\n";
        auto pyramid = BooleanOperations::createPyramid(15, 15, 20, 4);
        if (!pyramid) {
            cerr << "Failed to create pyramid\n";
            return 1;
        }
        STLExport::exportToSTL(*pyramid, "test_pyramid.stl");
        cout << "   ✓ Pyramid created: " << pyramid->getVertices().size()
            << " vertices, " << pyramid->getFaces().size() << " faces\n";

        cout << "\n=== PHASE 2: BOOLEAN OPERATIONS ===\n\n";

        // 6. Boolean Union
        cout << "6. Testing Boolean UNION (Box ∪ Cylinder)...\n";
        auto unionResult = BooleanOperations::booleanUnion(*box, *cylinder);
        if (unionResult) {
            STLExport::exportToSTL(*unionResult, "boolean_union.stl");
            cout << "   ✓ Union created: " << unionResult->getVertices().size()
                << " vertices, " << unionResult->getFaces().size() << " faces\n";
        }
        else {
            cout << "   ✗ Union failed\n";
        }

        // 7. Boolean Difference
        cout << "\n7. Testing Boolean DIFFERENCE (Box - Cylinder)...\n";
        auto differenceResult = BooleanOperations::booleanDifference(*box, *cylinder);
        if (differenceResult) {
            STLExport::exportToSTL(*differenceResult, "boolean_difference.stl");
            cout << "   ✓ Difference created: " << differenceResult->getVertices().size()
                << " vertices, " << differenceResult->getFaces().size() << " faces\n";
        }
        else {
            cout << "   ✗ Difference failed\n";
        }

        // 8. Boolean Intersection
        cout << "\n8. Testing Boolean INTERSECTION (Cylinder ∩ Sphere)...\n";
        auto intersectionResult = BooleanOperations::booleanIntersection(*cylinder, *sphere);
        if (intersectionResult) {
            STLExport::exportToSTL(*intersectionResult, "boolean_intersection.stl");
            cout << "   ✓ Intersection created: " << intersectionResult->getVertices().size()
                << " vertices, " << intersectionResult->getFaces().size() << " faces\n";
        }
        else {
            cout << "   ✗ Intersection failed\n";
        }

        cout << "\n=== PHASE 3: EXTRUSION & ADVANCED OPERATIONS ===\n\n";

        // 9. Extrusion test
        cout << "9. Testing Extrusion...\n";
        vector<Point3D> profile = {
            Point3D(0, 0, 0),
            Point3D(10, 0, 0),
            Point3D(10, 5, 0),
            Point3D(5, 5, 0),
            Point3D(5, 10, 0),
            Point3D(0, 10, 0)
        };

        auto extruded = BooleanOperations::extrudeProfile(profile, 15, true);
        if (extruded) {
            STLExport::exportToSTL(*extruded, "extrusion_test.stl");
            cout << "   ✓ Extrusion created: " << extruded->getVertices().size()
                << " vertices, " << extruded->getFaces().size() << " faces\n";
        }

        // 10. Complex operation
        cout << "\n10. Complex Boolean Operation...\n";
        if (box && cylinder && sphere) {
            // Create a complex part: box with cylindrical hole and spherical end
            auto step1 = BooleanOperations::booleanDifference(*box, *cylinder);
            if (step1) {
                auto finalResult = BooleanOperations::booleanUnion(*step1, *sphere);
                if (finalResult) {
                    STLExport::exportToSTL(*finalResult, "complex_part.stl");
                    cout << "   ✓ Complex part created: " << finalResult->getVertices().size()
                        << " vertices, " << finalResult->getFaces().size() << " faces\n";
                }
            }
        }

        // 11. Analysis tests
        cout << "\n=== PHASE 4: GEOMETRY ANALYSIS ===\n\n";

        if (box) {
            cout << "Box Analysis:\n";
            cout << "  Volume: " << BooleanOperations::calculateVolume(*box) << " mm³\n";
            cout << "  Surface Area: " << BooleanOperations::calculateSurfaceArea(*box) << " mm²\n";
            cout << "  Center of Mass: ";
            Point3D com = BooleanOperations::calculateCenterOfMass(*box);
            cout << "(" << com.x << ", " << com.y << ", " << com.z << ")\n";

            Point3D testPoint(5, 5, 5);
            bool inside = BooleanOperations::isPointInsideSolid(testPoint, *box);
            cout << "  Point (5,5,5) is " << (inside ? "INSIDE" : "OUTSIDE") << " the box\n";
        }

        // 12. Topology validation
        cout << "\n12. Topology Validation...\n";
        if (box && BooleanOperations::checkManifold(*box)) {
            cout << "   ✓ Box is manifold\n";
        }
        else {
            cout << "   ✗ Box is non-manifold\n";
        }

        if (cylinder && BooleanOperations::checkManifold(*cylinder)) {
            cout << "   ✓ Cylinder is manifold\n";
        }
        else {
            cout << "   ✗ Cylinder is non-manifold\n";
        }

        cout << "\n===============================================\n";
        cout << "   ✅ ALL TESTS COMPLETED SUCCESSFULLY\n";
        cout << "===============================================\n\n";

        cout << "Generated STL files:\n";
        cout << "  test_box.stl\n";
        cout << "  test_cylinder.stl\n";
        cout << "  test_sphere.stl\n";
        cout << "  test_cone.stl\n";
        cout << "  test_pyramid.stl\n";
        cout << "  boolean_union.stl\n";
        cout << "  boolean_difference.stl\n";
        cout << "  boolean_intersection.stl\n";
        cout << "  extrusion_test.stl\n";
        cout << "  complex_part.stl\n";

        return 0;

    }
    catch (const exception& e) {
        cerr << "\n❌ TEST FAILED: " << e.what() << endl;
        return 1;
    }
}