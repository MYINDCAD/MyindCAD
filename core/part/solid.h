// C:\MyIndustrialCAD\core\part\Solid.h
#pragma once

#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <algorithm>
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "Point3D.h"

class Solid {
private:
    std::string name_;
    std::vector<std::shared_ptr<Vertex>> vertices_;
    std::vector<std::shared_ptr<Edge>> edges_;
    std::vector<std::shared_ptr<Face>> faces_;
    Point3D bboxMin_;
    Point3D bboxMax_;

public:
    // Konstruktoren
    Solid() : name_("Unnamed Solid"), bboxMin_(0, 0, 0), bboxMax_(0, 0, 0) {}
    explicit Solid(const std::string& name) : name_(name), bboxMin_(0, 0, 0), bboxMax_(0, 0, 0) {}

    // Getter und Setter
    const std::string& getName() const { return name_; }
    void setName(const std::string& name) { name_ = name; }

    const std::vector<std::shared_ptr<Vertex>>& getVertices() const { return vertices_; }
    const std::vector<std::shared_ptr<Edge>>& getEdges() const { return edges_; }
    const std::vector<std::shared_ptr<Face>>& getFaces() const { return faces_; }

    // Elemente hinzufügen
    void addVertex(std::shared_ptr<Vertex> vertex) {
        vertices_.push_back(vertex);
        updateBoundingBox(vertex->getPosition());
    }

    void addEdge(std::shared_ptr<Edge> edge) {
        edges_.push_back(edge);
    }

    void addFace(std::shared_ptr<Face> face) {
        faces_.push_back(face);
    }

    // Elemente entfernen
    void clear() {
        vertices_.clear();
        edges_.clear();
        faces_.clear();
        bboxMin_ = bboxMax_ = Point3D(0, 0, 0);
    }

    // Bounding Box Management
    void updateBoundingBox(const Point3D& point) {
        if (vertices_.size() == 1) {
            // Erster Punkt
            bboxMin_ = bboxMax_ = point;
        }
        else {
            bboxMin_ = Point3D::min(bboxMin_, point);
            bboxMax_ = Point3D::max(bboxMax_, point);
        }
    }

    void getBoundingBox(Point3D& min, Point3D& max) const {
        min = bboxMin_;
        max = bboxMax_;
    }

    // Transformationen
    void translate(const Vector3D& translation) {
        for (auto& vertex : vertices_) {
            Point3D pos = vertex->getPosition();
            vertex->setPosition(pos + translation);
        }
        // Bounding Box aktualisieren
        bboxMin_ = bboxMin_ + translation;
        bboxMax_ = bboxMax_ + translation;
    }

    void scale(double factor) {
        for (auto& vertex : vertices_) {
            Point3D pos = vertex->getPosition();
            vertex->setPosition(Point3D(pos.x * factor, pos.y * factor, pos.z * factor));
        }
        // Bounding Box skalieren
        bboxMin_ = Point3D(bboxMin_.x * factor, bboxMin_.y * factor, bboxMin_.z * factor);
        bboxMax_ = Point3D(bboxMax_.x * factor, bboxMax_.y * factor, bboxMax_.z * factor);
    }

    // Volumenberechnung (vereinfacht)
    double volume() const {
        if (vertices_.empty()) return 0.0;

        // Einfache Berechnung für Quader
        double dx = bboxMax_.x - bboxMin_.x;
        double dy = bboxMax_.y - bboxMin_.y;
        double dz = bboxMax_.z - bboxMin_.z;

        return dx * dy * dz;
    }

    // Oberflächenberechnung (vereinfacht)
    double surfaceArea() const {
        if (faces_.empty()) return 0.0;

        double area = 0.0;
        // Hier würde echte Berechnung stehen
        return area;
    }

    // Validierung
    bool isValid() const {
        // Einfache Validierung
        return !vertices_.empty();
    }

    // Information
    void printInfo() const {
        std::cout << "Solid: " << name_ << std::endl;
        std::cout << "  Vertices: " << vertices_.size() << std::endl;
        std::cout << "  Edges: " << edges_.size() << std::endl;
        std::cout << "  Faces: " << faces_.size() << std::endl;
        std::cout << "  BBox Min: (" << bboxMin_.x << ", " << bboxMin_.y << ", " << bboxMin_.z << ")" << std::endl;
        std::cout << "  BBox Max: (" << bboxMax_.x << ", " << bboxMax_.y << ", " << bboxMax_.z << ")" << std::endl;
        std::cout << "  Volume: " << volume() << " mm³" << std::endl;
    }

    // Kopie erstellen
    std::shared_ptr<Solid> clone() const {
        auto clone = std::make_shared<Solid>(name_ + "_clone");

        // Vertices kopieren
        std::vector<std::shared_ptr<Vertex>> clonedVertices;
        for (const auto& v : vertices_) {
            auto clonedV = std::make_shared<Vertex>(v->getPosition());
            clonedVertices.push_back(clonedV);
            clone->addVertex(clonedV);
        }

        // Hier könnten Edges und Faces kopiert werden
        // (vereinfachte Implementierung)

        return clone;
    }
};