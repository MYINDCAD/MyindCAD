// C:\MyIndustrialCAD\core\part\Point3D.h
#pragma once
#include <algorithm>
#include <cmath>

struct Point3D {
    double x, y, z;

    Point3D(double x_ = 0, double y_ = 0, double z_ = 0)
        : x(x_), y(y_), z(z_) {
    }

    // Vergleichsoperator für std::map
    bool operator<(const Point3D& other) const {
        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        return z < other.z;
    }

    bool operator==(const Point3D& other) const {
        const double epsilon = 1e-10;
        return std::abs(x - other.x) < epsilon &&
            std::abs(y - other.y) < epsilon &&
            std::abs(z - other.z) < epsilon;
    }

    bool operator!=(const Point3D& other) const {
        return !(*this == other);
    }

    // Additionsoperator
    Point3D operator+(const Point3D& other) const {
        return Point3D(x + other.x, y + other.y, z + other.z);
    }

    // Subtraktionsoperator
    Point3D operator-(const Point3D& other) const {
        return Point3D(x - other.x, y - other.y, z - other.z);
    }

    // Zuweisungsoperatoren
    Point3D& operator=(const Point3D& other) {
        x = other.x;
        y = other.y;
        z = other.z;
        return *this;
    }

    // Min/Max Operatoren als statische Methoden
    static Point3D min(const Point3D& a, const Point3D& b) {
        return Point3D(
            std::min(a.x, b.x),
            std::min(a.y, b.y),
            std::min(a.z, b.z)
        );
    }

    static Point3D max(const Point3D& a, const Point3D& b) {
        return Point3D(
            std::max(a.x, b.x),
            std::max(a.y, b.y),
            std::max(a.z, b.z)
        );
    }
};

// Vector3D ist einfach ein Alias für Point3D
struct Vector3D : public Point3D {
    using Point3D::Point3D;
};