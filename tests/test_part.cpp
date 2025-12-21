// tests/test_part.cpp
#include <iostream>
#include "../core/part/Part.h"

int main() {
    std::cout << "=== MyIndustrialCAD Part Test ===\n" << std::endl;

    // 1. Part erstellen
    Part myPart("Test Part");
    std::cout << "1. Part erstellt: " << myPart.getName() << std::endl;

    // 2. Grundkörper erzeugen
    std::cout << "\n2. Grundkörper erzeugen..." << std::endl;

    auto box = Part::createBox(10, 20, 30);
    std::cout << "   - Box erstellt: " << box->getName() << std::endl;

    auto cylinder = Part::createCylinder(5, 15);
    std::cout << "   - Zylinder erstellt: " << cylinder->getName() << std::endl;

    // 3. Solids zum Part hinzufügen
    myPart.addSolid(box);
    myPart.addSolid(cylinder);
    std::cout << "\n3. Solids hinzugefügt. Anzahl: " << myPart.getSolidCount() << std::endl;

    // 4. Bounding Box abfragen
    Point3D min, max;
    myPart.getBoundingBox(min, max);
    std::cout << "\n4. Bounding Box:" << std::endl;
    std::cout << "   Min: (" << min.x << ", " << min.y << ", " << min.z << ")" << std::endl;
    std::cout << "   Max: (" << max.x << ", " << max.y << ", " << max.z << ")" << std::endl;

    // 5. Volumen berechnen
    double volume = myPart.totalVolume();
    std::cout << "\n5. Gesamtvolumen: " << volume << " mm³" << std::endl;

    // 6. Validierung
    std::cout << "\n6. Validierung:" << std::endl;
    std::cout << "   Part ist " << (myPart.isValid() ? "gültig" : "ungültig") << std::endl;

    // 7. Boolean Operations testen
    std::cout << "\n7. Boolean Operations testen..." << std::endl;
    auto unionResult = Part::booleanUnion(*box, *cylinder);
    std::cout << "   Union durchgeführt" << std::endl;

    auto diffResult = Part::booleanDifference(*box, *cylinder);
    std::cout << "   Difference durchgeführt" << std::endl;

    auto intersectResult = Part::booleanIntersection(*box, *cylinder);
    std::cout << "   Intersection durchgeführt" << std::endl;

    std::cout << "\n=== Test erfolgreich abgeschlossen ===" << std::endl;

    return 0;
}