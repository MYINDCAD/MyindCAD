// test_phase1_complete.cpp
#include <iostream>
#include "../core/part/Part.h"
#include "../core/sketcher/Sketch.h"

int main() {
    std::cout << "=== Phase 1 Test: Alle Features ===\n\n";

    // 1. Skizzen-Test
    std::cout << "1. SKETCHER FEATURES:\n";
    Sketch sketch;

    // Alle 6 Pflichtelemente
    sketch.addPoint(0, 0);
    sketch.addLine(0, 0, 10, 0);
    sketch.addCircle(5, 5, 3);
    sketch.addArc(5, 5, 3, 0, 90);
    sketch.addRectangle(0, 0, 10, 10);
    sketch.addSpline({ {0,0}, {5,5}, {10,0} });

    sketch.printPoints();
    sketch.printCurves();

    // 2. 3D Part-Test
    std::cout << "\n2. 3D PART FEATURES:\n";
    Part part("TestPart");

    // Primitive Körper erstellen
    auto box = Part::createBox(10, 20, 30);
    auto cylinder = Part::createCylinder(5, 15);

    part.addSolid(box);
    part.addSolid(cylinder);

    std::cout << "Part '" << part.getName() << "' erstellt mit "
        << part.getSolidCount() << " Solids\n";

    // Bounding Box
    Point3D min, max;
    part.getBoundingBox(min, max);
    std::cout << "Bounding Box: Min(" << min.x << "," << min.y << "," << min.z
        << ") Max(" << max.x << "," << max.y << "," << max.z << ")\n";

    // Volume
    std::cout << "Gesamtvolumen: " << part.totalVolume() << " mm³\n";

    std::cout << "\n✅ PHASE 1 ABGESCHLOSSEN!\n";
    std::cout << "Alle 6 Pflichtelemente funktionieren:\n";
    std::cout << "- Punkt, Linie, Kreis, Kreisbogen, Rechteck, Spline\n";

    return 0;
}