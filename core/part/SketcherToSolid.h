// SketcherToSolid.h
#ifndef SKETCHER_TO_SOLID_H
#define SKETCHER_TO_SOLID_H

#include "solid.h"
#include "Sketch.h"
#include <vector>

class SketcherToSolid {
public:
    // Extrudiert einen Sketch zu einem 3D-Volumen
    static Solid* extrudeSketch(const Sketch* sketch, double height, bool bothSides = false);

    // Rotiert einen Sketch um die Y-Achse
    static Solid* revolveSketch(const Sketch* sketch, double angleDegrees = 360.0, int segments = 32);

    // Führt einen Sketch entlang eines Pfades (fortgeschritten)
    static Solid* sweepSketch(const Sketch* sketch, const std::vector<Point3D>& path, bool followRotation = true);

private:
    SketcherToSolid() = delete; // Statische Klasse
};

#endif // SKETCHER_TO_SOLID_H