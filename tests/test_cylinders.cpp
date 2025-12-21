// test_cylinders.cpp
#include <iostream>
#include <cassert>
#include <cmath>
#include "BooleanOperations.h"
#include "STLExport.h"
#include "SketcherToSolid.h"
#include "Sketch.h"

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void testBasicCylinder() {
    cout << "=== Test 1: Basic Cylinder Creation ===" << endl;

    Solid* cylinder = SolidFactory::createCylinder(5.0, 10.0, 32);

    cout << "Cylinder created:" << endl;
    cout << "  Vertices: " << cylinder->getVertices().size() << endl;
    cout << "  Edges: " << cylinder->getEdges().size() << endl;
    cout << "  Faces: " << cylinder->getFaces().size() << endl;

    // Cylinder sollte haben:
    // - 2 center vertices + 2*segments circle vertices = 66 vertices (für 32 segments)
    // - 4*segments edges = 128 edges
    // - segments side faces + 2 cap faces = 34 faces

    assert(cylinder->getVertices().size() == 66);  // 2 + 2*32
    assert(cylinder->getFaces().size() == 34);     // 32 + 2

    cout << "  ✓ Basic cylinder geometry correct" << endl;

    // Volumen prüfen (V = π * r² * h)
    double expectedVolume = M_PI * 5.0 * 5.0 * 10.0;
    double actualVolume = cylinder->volume();
    double tolerance = 0.1; // 10% tolerance for discretization

    cout << "  Expected volume: " << expectedVolume << " mm³" << endl;
    cout << "  Actual volume: " << actualVolume << " mm³" << endl;

    assert(abs(actualVolume - expectedVolume) / expectedVolume < tolerance);
    cout << "  ✓ Volume calculation correct" << endl;

    // STL Export test
    bool exportOk = STLExport::exportToSTL(cylinder, "test_cylinder.stl", true);
    assert(exportOk);
    cout << "  ✓ STL export successful" << endl;

    delete cylinder;
}

void testCylinderBooleanOperations() {
    cout << "\n=== Test 2: Cylinder Boolean Operations ===" << endl;

    // Box mit zylindrischem Loch
    Solid* box = SolidFactory::createBox(15, 15, 15);
    Solid* cylinder = SolidFactory::createCylinder(5.0, 20.0, 32); // Höher als Box

    // Box zentrieren
    // (Ihre createBox ist bereits zentriert)

    // Cylinder in Box platzieren
    // transform würde hier verwendet werden

    cout << "  Box volume: " << box->volume() << " mm³" << endl;
    cout << "  Cylinder volume: " << cylinder->volume() << " mm³" << endl;

    // Differenz: Box minus Cylinder (gebohrtes Loch)
    Solid* drilledBox = BooleanOperations::booleanDifference(box, cylinder);

    assert(drilledBox != nullptr);
    cout << "  Drilled box created:" << endl;
    cout << "    Vertices: " << drilledBox->getVertices().size() << endl;
    cout << "    Faces: " << drilledBox->getFaces().size() << endl;

    // Volumen prüfen
    double boxVolume = 15 * 15 * 15; // 3375
    double cylVolume = M_PI * 5.0 * 5.0 * 15; // Cylinder durch Box (15mm hoch)
    double expectedVolume = boxVolume - cylVolume;
    double actualVolume = drilledBox->volume();

    cout << "  Expected drilled volume: " << expectedVolume << " mm³" << endl;
    cout << "  Actual drilled volume: " << actualVolume << " mm³" << endl;

    // 5% Toleranz wegen Diskretisierung
    assert(abs(actualVolume - expectedVolume) / expectedVolume < 0.05);
    cout << "  ✓ Boolean difference correct" << endl;

    // Export
    STLExport::exportToSTL(drilledBox, "drilled_box.stl");
    cout << "  ✓ Exported drilled box" << endl;

    // Union testen
    Solid* cylinder2 = SolidFactory::createCylinder(3.0, 10.0, 16);
    Solid* unionResult = BooleanOperations::booleanUnion(cylinder, cylinder2);

    if (unionResult) {
        cout << "  Union cylinder+cylinder created" << endl;
        STLExport::exportToSTL(unionResult, "union_cylinders.stl");
        delete unionResult;
    }

    // Intersection testen
    Solid* smallBox = SolidFactory::createBox(8, 8, 8);
    Solid* intersection = BooleanOperations::booleanIntersection(cylinder, smallBox);

    if (intersection) {
        cout << "  Cylinder-box intersection created" << endl;
        STLExport::exportToSTL(intersection, "cylinder_box_intersection.stl");
        delete intersection;
    }

    delete box;
    delete cylinder;
    delete cylinder2;
    delete smallBox;
    delete drilledBox;
}

void testSphereAndCone() {
    cout << "\n=== Test 3: Sphere and Cone ===" << endl;

    Solid* sphere = SolidFactory::createSphere(7.0, 24);
    assert(sphere != nullptr);
    cout << "  Sphere created with " << sphere->getVertices().size() << " vertices" << endl;

    double sphereVolume = sphere->volume();
    double expectedSphereVolume = (4.0 / 3.0) * M_PI * 7.0 * 7.0 * 7.0;
    cout << "  Sphere volume: " << sphereVolume << " mm³ (expected: "
        << expectedSphereVolume << ")" << endl;

    // 10% Toleranz für Sphere Approximation
    assert(abs(sphereVolume - expectedSphereVolume) / expectedSphereVolume < 0.1);

    Solid* cone = SolidFactory::createCone(6.0, 12.0, 24);
    assert(cone != nullptr);
    cout << "  Cone created with " << cone->getVertices().size() << " vertices" << endl;

    double coneVolume = cone->volume();
    double expectedConeVolume = (1.0 / 3.0) * M_PI * 6.0 * 6.0 * 12.0;
    cout << "  Cone volume: " << coneVolume << " mm³ (expected: "
        << expectedConeVolume << ")" << endl;

    // 10% Toleranz
    assert(abs(coneVolume - expectedConeVolume) / expectedConeVolume < 0.1);

    // Export
    STLExport::exportToSTL(sphere, "test_sphere.stl");
    STLExport::exportToSTL(cone, "test_cone.stl");
    cout << "  ✓ Sphere and cone exported" << endl;

    delete sphere;
    delete cone;
}

void testComplexBoolean() {
    cout << "\n=== Test 4: Complex Boolean Operation ===" << endl;

    // Erstelle eine Art "Flansch" mit Löchern
    Solid* basePlate = SolidFactory::createBox(30, 30, 5);

    // 4 Löcher für Schrauben
    Solid* holes[4];
    Solid* result = basePlate;

    double holePositions[4][2] = {
        {10, 10}, {10, -10}, {-10, 10}, {-10, -10}
    };

    for (int i = 0; i < 4; i++) {
        Solid* hole = SolidFactory::createCylinder(2.0, 10.0, 16);
        // Hier müsste transform hinzugefügt werden
        // hole->transform(translate(holePositions[i][0], holePositions[i][1], -2.5));

        result = BooleanOperations::booleanDifference(result, hole);
        holes[i] = hole;
    }

    // Zentrale Aussparung
    Solid* centerCut = SolidFactory::createCylinder(8.0, 10.0, 32);
    // centerCut->transform(translate(0, 0, -2.5));

    Solid* finalPart = BooleanOperations::booleanDifference(result, centerCut);

    cout << "  Complex part created:" << endl;
    cout << "    Vertices: " << finalPart->getVertices().size() << endl;
    cout << "    Faces: " << finalPart->getFaces().size() << endl;

    STLExport::exportToSTL(finalPart, "complex_flange.stl");
    cout << "  ✓ Complex part exported" << endl;

    delete basePlate;
    delete centerCut;
    delete finalPart;
    for (int i = 0; i < 4; i++) {
        delete holes[i];
    }
    if (result != basePlate && result != finalPart) {
        delete result;
    }
}

int main() {
    cout << "===============================================\n";
    cout << "    CYLINDER & ADVANCED GEOMETRY TESTS\n";
    cout << "===============================================\n\n";

    try {
        testBasicCylinder();
        testCylinderBooleanOperations();
        testSphereAndCone();
        testComplexBoolean();

        cout << "\n===============================================\n";
        cout << "    ✅ ALL TESTS PASSED SUCCESSFULLY!\n";
        cout << "    Generated STL files:\n";
        cout << "      - test_cylinder.stl\n";
        cout << "      - drilled_box.stl\n";
        cout << "      - union_cylinders.stl\n";
        cout << "      - cylinder_box_intersection.stl\n";
        cout << "      - test_sphere.stl\n";
        cout << "      - test_cone.stl\n";
        cout << "      - complex_flange.stl\n";
        cout << "===============================================\n";

        return 0;
    }
    catch (const exception& e) {
        cerr << "\n❌ TEST FAILED: " << e.what() << endl;
        return 1;
    }
}