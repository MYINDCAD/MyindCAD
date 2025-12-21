// C:\MyIndustrialCAD\core\part\Part.cpp
#include "Part.h"
#include "BooleanOperations.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <memory>
#include <iostream>
#include <cmath>

// Sketch.h hier includieren (nicht in Part.h)
#include "../sketcher/Sketch.h"

// PI für Windows definieren
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {
    std::shared_ptr<Vertex> createVertex(double x, double y, double z) {
        return std::make_shared<Vertex>(Point3D(x, y, z));
    }

    std::shared_ptr<Edge> createEdge(std::shared_ptr<Vertex> v1,
        std::shared_ptr<Vertex> v2) {
        return std::make_shared<Edge>(v1, v2);
    }

    std::shared_ptr<Face> createFace(const std::vector<std::shared_ptr<Edge>>& edges,
        const std::vector<std::shared_ptr<Vertex>>& vertices) {
        return std::make_shared<Face>(edges, vertices);
    }
}

// Konstruktoren
Part::Part() : name_("Unnamed Part") {}

Part::Part(const std::string& name) : name_(name) {}

// Solid Management
void Part::addSolid(std::shared_ptr<Solid> solid) {
    solids_.push_back(solid);
}

std::shared_ptr<Solid> Part::getSolid(size_t index) const {
    if (index < solids_.size()) {
        return solids_[index];
    }
    return nullptr;
}

void Part::removeSolid(int index) {
    if (index >= 0 && index < static_cast<int>(solids_.size())) {
        solids_.erase(solids_.begin() + index);
    }
}

// Primitive Solids
std::shared_ptr<Solid> Part::createBox(double width, double height, double depth) {
    auto solid = std::make_shared<Solid>("Box");

    // 8 Vertices
    std::vector<std::shared_ptr<Vertex>> vertices = {
        createVertex(0, 0, 0),
        createVertex(width, 0, 0),
        createVertex(width, height, 0),
        createVertex(0, height, 0),
        createVertex(0, 0, depth),
        createVertex(width, 0, depth),
        createVertex(width, height, depth),
        createVertex(0, height, depth)
    };

    // 12 Edges
    std::vector<std::shared_ptr<Edge>> edges = {
        createEdge(vertices[0], vertices[1]),
        createEdge(vertices[1], vertices[2]),
        createEdge(vertices[2], vertices[3]),
        createEdge(vertices[3], vertices[0]),
        createEdge(vertices[4], vertices[5]),
        createEdge(vertices[5], vertices[6]),
        createEdge(vertices[6], vertices[7]),
        createEdge(vertices[7], vertices[4]),
        createEdge(vertices[0], vertices[4]),
        createEdge(vertices[1], vertices[5]),
        createEdge(vertices[2], vertices[6]),
        createEdge(vertices[3], vertices[7])
    };

    // 6 Faces
    std::vector<std::shared_ptr<Face>> faces = {
        createFace({edges[0], edges[1], edges[2], edges[3]},
                   {vertices[0], vertices[1], vertices[2], vertices[3]}),
        createFace({edges[4], edges[5], edges[6], edges[7]},
                   {vertices[4], vertices[5], vertices[6], vertices[7]}),
        createFace({edges[0], edges[9], edges[4], edges[8]},
                   {vertices[0], vertices[1], vertices[5], vertices[4]}),
        createFace({edges[2], edges[10], edges[6], edges[11]},
                   {vertices[2], vertices[3], vertices[7], vertices[6]}),
        createFace({edges[3], edges[11], edges[7], edges[8]},
                   {vertices[3], vertices[0], vertices[4], vertices[7]}),
        createFace({edges[1], edges[10], edges[5], edges[9]},
                   {vertices[1], vertices[2], vertices[6], vertices[5]})
    };

    for (auto& v : vertices) solid->addVertex(v);
    for (auto& e : edges) solid->addEdge(e);
    for (auto& f : faces) solid->addFace(f);

    return solid;
}

std::shared_ptr<Solid> Part::createCylinder(double radius, double height, int segments) {
    auto solid = std::make_shared<Solid>("Cylinder");
    if (segments < 3) segments = 3;

    // Vertices
    auto bottomCenter = createVertex(0, 0, 0);
    auto topCenter = createVertex(0, 0, height);
    solid->addVertex(bottomCenter);
    solid->addVertex(topCenter);

    std::vector<std::shared_ptr<Vertex>> bottomVertices;
    std::vector<std::shared_ptr<Vertex>> topVertices;

    for (int i = 0; i < segments; ++i) {
        double angle = 2.0 * M_PI * i / segments;
        double x = radius * std::cos(angle);
        double y = radius * std::sin(angle);

        auto bottomVert = createVertex(x, y, 0);
        auto topVert = createVertex(x, y, height);

        bottomVertices.push_back(bottomVert);
        topVertices.push_back(topVert);
        solid->addVertex(bottomVert);
        solid->addVertex(topVert);
    }

    std::cout << "Cylinder created (basic structure)" << std::endl;
    return solid;
}

std::shared_ptr<Solid> Part::createSphere(double radius, int segments) {
    auto solid = std::make_shared<Solid>("Sphere");
    std::cout << "Sphere placeholder created" << std::endl;
    return solid;
}

std::shared_ptr<Solid> Part::createCone(double bottomRadius, double topRadius, double height, int segments) {
    auto solid = std::make_shared<Solid>("Cone");
    std::cout << "Cone placeholder created" << std::endl;
    return solid;
}

// Sketch-basierte Solids - Implementierungen
std::shared_ptr<Solid> Part::extrudeSketch(const Sketch& sketch, double height, bool bothSides) {
    std::cout << "Extruding sketch to height: " << height << " mm" << std::endl;

    auto solid = std::make_shared<Solid>("Extruded");

    // Platzhalter-Implementierung
    // In der realen Implementierung würden Sie die Sketch-Punkte verwenden

    // Einfache Box als Platzhalter erstellen
    auto placeholderBox = createBox(10, 10, height);

    // Eigenschaften kopieren (vereinfacht)
    return placeholderBox;
}

std::shared_ptr<Solid> Part::revolveSketch(const Sketch& sketch, double angleDegrees,
    const Point3D& axisPoint, const Vector3D& axisDir) {
    std::cout << "Revolving sketch by " << angleDegrees << " degrees" << std::endl;

    auto solid = std::make_shared<Solid>("Revolved");

    // Platzhalter-Implementierung
    auto placeholderCylinder = createCylinder(5, 10);

    return placeholderCylinder;
}

std::shared_ptr<Solid> Part::createSolidFromSketchProfile(const Sketch& sketch, double z) {
    std::cout << "Creating solid from sketch profile at z=" << z << std::endl;

    auto solid = std::make_shared<Solid>("SketchProfile");

    // Platzhalter
    auto placeholderBox = createBox(5, 5, 1);

    return placeholderBox;
}

// Hilfsfunktionen
std::vector<Point3D> Part::convert2DPointsTo3D(const std::vector<SimplePoint>& points2D, double z) {
    std::vector<Point3D> points3D;
    for (const auto& p : points2D) {
        points3D.emplace_back(p.x, p.y, z);
    }
    return points3D;
}

// Boolean Operations
std::shared_ptr<Solid> Part::booleanUnion(const Solid& a, const Solid& b) {
    return BooleanOperations::booleanUnion(a, b);
}

std::shared_ptr<Solid> Part::booleanDifference(const Solid& a, const Solid& b) {
    return BooleanOperations::booleanDifference(a, b);
}

std::shared_ptr<Solid> Part::booleanIntersection(const Solid& a, const Solid& b) {
    return BooleanOperations::booleanIntersection(a, b);
}

// Transformationen
void Part::translate(double x, double y, double z) {
    Vector3D translation(x, y, z);
    for (auto& solid : solids_) {
        solid->translate(translation);
    }
}

void Part::rotate(double angleDegrees, double axisX, double axisY, double axisZ) {
    std::cout << "Rotating part by " << angleDegrees << " degrees" << std::endl;
    // Platzhalter
}

void Part::scale(double factor) {
    std::cout << "Scaling part by factor " << factor << std::endl;
    // Platzhalter
}

// Abfragen
double Part::totalVolume() const {
    double volume = 0.0;
    for (const auto& solid : solids_) {
        volume += solid->volume();
    }
    return volume;
}

void Part::getBoundingBox(Point3D& min, Point3D& max) const {
    if (solids_.empty()) {
        min = max = Point3D(0, 0, 0);
        return;
    }

    solids_[0]->getBoundingBox(min, max);

    for (size_t i = 1; i < solids_.size(); ++i) {
        Point3D solidMin, solidMax;
        solids_[i]->getBoundingBox(solidMin, solidMax);

        min.x = std::min(min.x, solidMin.x);
        min.y = std::min(min.y, solidMin.y);
        min.z = std::min(min.z, solidMin.z);

        max.x = std::max(max.x, solidMax.x);
        max.y = std::max(max.y, solidMax.y);
        max.z = std::max(max.z, solidMax.z);
    }
}

// Export
bool Part::exportToSTEP(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) return false;

    file << "ISO-10303-21;\n";
    file << "HEADER;\n";
    file << "FILE_DESCRIPTION(('MyIndustrialCAD Part'),'2;1');\n";
    file << "FILE_NAME('" << name_ << "','2024-12-15T12:00:00','','','MyIndustrialCAD','','');\n";
    file << "FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));\n";
    file << "ENDSEC;\n";
    file << "DATA;\n";
    file << "ENDSEC;\n";
    file << "END-ISO-10303-21;\n";

    file.close();
    std::cout << "Exported STEP: " << filename << std::endl;
    return true;
}

bool Part::exportToSTL(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) return false;

    file << "solid " << name_ << "\n";
    file << "endsolid " << name_ << "\n";

    file.close();
    std::cout << "Exported STL: " << filename << std::endl;
    return true;
}

bool Part::exportToOBJ(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) return false;

    file << "# MyIndustrialCAD OBJ Export\n";
    file << "# Part: " << name_ << "\n";

    file.close();
    std::cout << "Exported OBJ: " << filename << std::endl;
    return true;
}

// Validierung
bool Part::isValid() const {
    for (const auto& solid : solids_) {
        if (solid->getVertices().empty()) return false;
    }
    return true;
}

void Part::repair() {
    std::cout << "Part repair function called" << std::endl;
}