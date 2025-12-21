// C:\MyIndustrialCAD\core\part\Vertex.h
#pragma once

#include "Point3D.h"  // WICHTIG: Diese Datei muss existieren

class Vertex {
private:
    Point3D position_;  // Hier war der Fehler: Typspezifizierer fehlte

public:
    // Konstruktor
    Vertex(const Point3D& pos = Point3D(0, 0, 0)) : position_(pos) {}

    // Getter
    Point3D getPosition() const { return position_; }

    // Setter
    void setPosition(const Point3D& pos) { position_ = pos; }

    // Transformation
    void translate(const Vector3D& v) {
        position_ = position_ + v;
    }

    // Vergleichsoperator
    bool operator==(const Vertex& other) const {
        return position_ == other.position_;
    }
};