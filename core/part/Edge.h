// C:\MyIndustrialCAD\core\part\Edge.h
#pragma once
#include <memory>
#include <cmath>  // HIER HINZUFÜGEN!
#include "Vertex.h"

class Edge {
private:
    std::shared_ptr<Vertex> v1_, v2_;

public:
    Edge(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2)
        : v1_(v1), v2_(v2) {
    }

    std::shared_ptr<Vertex> getV1() const { return v1_; }
    std::shared_ptr<Vertex> getV2() const { return v2_; }

    double length() const {
        auto p1 = v1_->getPosition();
        auto p2 = v2_->getPosition();
        double dx = p2.x - p1.x;
        double dy = p2.y - p1.y;
        double dz = p2.z - p1.z;
        return std::sqrt(dx * dx + dy * dy + dz * dz);  // std::sqrt statt sqrt
    }
};