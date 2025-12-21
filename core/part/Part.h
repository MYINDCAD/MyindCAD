// C:\MyIndustrialCAD\core\part\Part.h
#pragma once

#include <memory>
#include <vector>
#include <string>
#include "Solid.h"

// Forward declarations für Sketch und Point
class Sketch;
struct Point;

// Simple Point struct für 2D (falls nicht in Sketch definiert)
struct SimplePoint {
    double x, y;
    SimplePoint(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}
};

class Part {
private:
    std::string name_;
    std::vector<std::shared_ptr<Solid>> solids_;

public:
    // Konstruktoren
    Part();
    explicit Part(const std::string& name);

    // Getter und Setter
    const std::string& getName() const { return name_; }
    void setName(const std::string& name) { name_ = name; }

    // Solid Management
    void addSolid(std::shared_ptr<Solid> solid);
    std::shared_ptr<Solid> getSolid(size_t index) const;
    size_t getSolidCount() const { return solids_.size(); }
    void removeSolid(int index);

    const std::vector<std::shared_ptr<Solid>>& getSolids() const { return solids_; }

    // Primitive Solids
    static std::shared_ptr<Solid> createBox(double width, double height, double depth);
    static std::shared_ptr<Solid> createCylinder(double radius, double height, int segments = 32);
    static std::shared_ptr<Solid> createSphere(double radius, int segments = 16);
    static std::shared_ptr<Solid> createCone(double bottomRadius, double topRadius,
        double height, int segments = 32);

    // Sketch-basierte Solids (mit Forward Declared Sketch)
    static std::shared_ptr<Solid> extrudeSketch(const Sketch& sketch, double height,
        bool bothSides = false);

    static std::shared_ptr<Solid> revolveSketch(const Sketch& sketch, double angleDegrees,
        const Point3D& axisPoint = Point3D(0, 0, 0),
        const Vector3D& axisDir = Vector3D(0, 0, 1));

    // Boolean Operations
    static std::shared_ptr<Solid> booleanUnion(const Solid& a, const Solid& b);
    static std::shared_ptr<Solid> booleanDifference(const Solid& a, const Solid& b);
    static std::shared_ptr<Solid> booleanIntersection(const Solid& a, const Solid& b);

    // Transformationen
    void translate(double x, double y, double z);
    void rotate(double angleDegrees, double axisX, double axisY, double axisZ);
    void scale(double factor);

    // Hilfsfunktionen
    static std::vector<Point3D> convert2DPointsTo3D(const std::vector<SimplePoint>& points2D, double z = 0);

    static std::shared_ptr<Solid> createSolidFromSketchProfile(const Sketch& sketch, double z = 0);

    // Abfragen
    double totalVolume() const;
    void getBoundingBox(Point3D& min, Point3D& max) const;

    // Export
    bool exportToSTEP(const std::string& filename) const;
    bool exportToSTL(const std::string& filename) const;
    bool exportToOBJ(const std::string& filename) const;

    // Validierung
    bool isValid() const;
    void repair();
};