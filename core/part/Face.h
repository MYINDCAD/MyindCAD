// C:\MyIndustrialCAD\core\part\Face.h
#pragma once
#include <vector>
#include <memory>
#include "Edge.h"
#include "Vertex.h"

class Face {
private:
    std::vector<std::shared_ptr<Edge>> edges_;
    std::vector<std::shared_ptr<Vertex>> vertices_;

public:
    Face(const std::vector<std::shared_ptr<Edge>>& edges,
        const std::vector<std::shared_ptr<Vertex>>& vertices)
        : edges_(edges), vertices_(vertices) {
    }

    const std::vector<std::shared_ptr<Edge>>& getEdges() const { return edges_; }
    const std::vector<std::shared_ptr<Vertex>>& getVertices() const { return vertices_; }

    bool isPlanar() const {
        // Vereinfachte Implementierung
        return vertices_.size() >= 3;
    }
};