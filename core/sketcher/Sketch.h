// C:\MyIndustrialCAD\sketcher\Sketch.h
#pragma once
#include <vector>

struct Point {
    double x, y;
    Point(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}
};

class Sketch {
public:
    Sketch() = default;

    std::vector<Point> getBoundaryPoints() const {
        return std::vector<Point>{Point(0, 0), Point(10, 0), Point(10, 10), Point(0, 10)};
    }

    bool isValidForExtrusion() const { return true; }
};