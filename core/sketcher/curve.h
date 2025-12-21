#pragma once
#include "Point.h"

enum class CurveType {
    Line,
    Circle,
    Arc,
    Rectangle,
    Spline
};

struct Curve {
    CurveType type;

    // Linien / Rechtecke
    Point start;
    Point end;

    // Kreise / Bögen
    Point center;
    double radius = 0.0;
    double startAngle = 0.0; // für Arc
    double endAngle = 0.0;   // für Arc

    // Rechteck-spezifisch
    Point corner2;

    // Spline (ein Typ)
    std::vector<Point> controlPoints;
};
