#include "Sketch.h"
#include <iostream>
#include <algorithm>
#include <cmath>

// PI Konstante definieren (für Windows, da M_PI nicht standard ist)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// === Punkte ===
void Sketch::addPoint(double x, double y) {
    points_.push_back({ x, y });
}

void Sketch::removePoint(int index) {
    if (validPointIndex(index)) {
        points_.erase(points_.begin() + index);
    }
}

void Sketch::clearPoints() {
    points_.clear();
}

const std::vector<Point>& Sketch::getPoints() const {
    return points_;
}

void Sketch::printPoints() const {
    std::cout << "Points (" << points_.size() << "):\n";
    for (size_t i = 0; i < points_.size(); ++i)
        std::cout << "  " << i << ": (" << points_[i].x << ", " << points_[i].y << ")\n";
}

// === Kurven ===
void Sketch::addLine(double x1, double y1, double x2, double y2) {
    Curve c;
    c.type = CurveType::Line;
    c.start = { x1, y1 };
    c.end = { x2, y2 };
    curves_.push_back(c);
}

void Sketch::addCircle(double centerX, double centerY, double radius) {
    Curve c;
    c.type = CurveType::Circle;
    c.center = { centerX, centerY };
    c.radius = radius;
    curves_.push_back(c);
}

void Sketch::addArc(double centerX, double centerY, double radius,
    double startAngle, double endAngle) {
    Curve c;
    c.type = CurveType::Arc;
    c.center = { centerX, centerY };
    c.radius = radius;
    c.startAngle = startAngle;
    c.endAngle = endAngle;
    curves_.push_back(c);
}

void Sketch::addRectangle(double x1, double y1, double x2, double y2) {
    Curve c;
    c.type = CurveType::Rectangle;
    c.start = { x1, y1 };
    c.corner2 = { x2, y2 };
    curves_.push_back(c);
}

void Sketch::addSpline(const std::vector<Point>& controlPoints) {
    if (controlPoints.empty()) return;
    Curve c;
    c.type = CurveType::Spline;
    c.controlPoints = controlPoints;
    curves_.push_back(c);
}

void Sketch::removeCurve(int index) {
    if (validCurveIndex(index)) {
        curves_.erase(curves_.begin() + index);
    }
}

void Sketch::clearCurves() {
    curves_.clear();
}

const std::vector<Curve>& Sketch::getCurves() const {
    return curves_;
}

void Sketch::printCurves() const {
    std::cout << "Curves (" << curves_.size() << "):\n";
    for (size_t i = 0; i < curves_.size(); ++i) {
        const Curve& c = curves_[i];
        std::cout << "  " << i << ": ";
        switch (c.type) {
        case CurveType::Line:
            std::cout << "Line (" << c.start.x << "," << c.start.y
                << ") -> (" << c.end.x << "," << c.end.y << ")";
            break;
        case CurveType::Circle:
            std::cout << "Circle center(" << c.center.x << "," << c.center.y
                << "), r=" << c.radius;
            break;
        case CurveType::Arc:
            std::cout << "Arc center(" << c.center.x << "," << c.center.y
                << "), r=" << c.radius << ", start=" << c.startAngle
                << "°, end=" << c.endAngle << "°";
            break;
        case CurveType::Rectangle:
            std::cout << "Rectangle (" << c.start.x << "," << c.start.y
                << ") -> (" << c.corner2.x << "," << c.corner2.y << ")";
            break;
        case CurveType::Spline:
            std::cout << "Spline with " << c.controlPoints.size() << " control points";
            break;
        }
        std::cout << "\n";
    }
}

// === Constraints ===
bool Sketch::validCurveIndex(int index) const {
    return index >= 0 && index < static_cast<int>(curves_.size());
}

bool Sketch::validPointIndex(int index) const {
    return index >= 0 && index < static_cast<int>(points_.size());
}

bool Sketch::canAddConstraint(const Constraint& c) const {
    if (c.entityA < 0) return false;
    if (c.entityB >= 0 && c.entityA == c.entityB) return false;

    // Prüfe ob Entities existieren
    bool entityAValid = (c.entityA < static_cast<int>(curves_.size()) ||
        c.entityA < static_cast<int>(points_.size()));
    bool entityBValid = (c.entityB < 0 ||
        c.entityB < static_cast<int>(curves_.size()) ||
        c.entityB < static_cast<int>(points_.size()));

    return entityAValid && entityBValid;
}

bool Sketch::addHorizontalConstraint(int curveIndex) {
    if (!validCurveIndex(curveIndex)) return false;
    constraints_.push_back({ ConstraintType::Horizontal, curveIndex, -1, 0.0 });
    return true;
}

bool Sketch::addVerticalConstraint(int curveIndex) {
    if (!validCurveIndex(curveIndex)) return false;
    constraints_.push_back({ ConstraintType::Vertical, curveIndex, -1, 0.0 });
    return true;
}

bool Sketch::addDistanceConstraint(int entityA, int entityB, double distance) {
    if ((!validCurveIndex(entityA) && !validPointIndex(entityA)) ||
        (!validCurveIndex(entityB) && !validPointIndex(entityB))) {
        return false;
    }
    constraints_.push_back({ ConstraintType::Distance, entityA, entityB, distance });
    return true;
}

bool Sketch::addAngleConstraint(int curveA, int curveB, double angleDegrees) {
    if (!validCurveIndex(curveA) || !validCurveIndex(curveB)) return false;
    constraints_.push_back({ ConstraintType::Angle, curveA, curveB, angleDegrees });
    return true;
}

bool Sketch::addTangentConstraint(int curveA, int curveB) {
    if (!validCurveIndex(curveA) || !validCurveIndex(curveB)) return false;
    constraints_.push_back({ ConstraintType::Tangent, curveA, curveB, 0.0 });
    return true;
}

bool Sketch::addCoincidentConstraint(int pointIndex, int curveIndex) {
    if (!validPointIndex(pointIndex) || !validCurveIndex(curveIndex)) return false;
    constraints_.push_back({ ConstraintType::Coincident, pointIndex, curveIndex, 0.0 });
    return true;
}

void Sketch::removeConstraint(int index) {
    if (index >= 0 && index < static_cast<int>(constraints_.size())) {
        constraints_.erase(constraints_.begin() + index);
    }
}

void Sketch::clearConstraints() {
    constraints_.clear();
}

const std::vector<Constraint>& Sketch::getConstraints() const {
    return constraints_;
}

void Sketch::printConstraints() const {
    std::cout << "Constraints (" << constraints_.size() << "):\n";
    for (size_t i = 0; i < constraints_.size(); ++i) {
        const Constraint& c = constraints_[i];
        std::cout << "  " << i << ": ";
        switch (c.type) {
        case ConstraintType::Horizontal:
            std::cout << "Horizontal on curve " << c.entityA;
            break;
        case ConstraintType::Vertical:
            std::cout << "Vertical on curve " << c.entityA;
            break;
        case ConstraintType::Distance:
            std::cout << "Distance between " << c.entityA << " and "
                << c.entityB << " = " << c.value;
            break;
        case ConstraintType::Angle:
            std::cout << "Angle between " << c.entityA << " and "
                << c.entityB << " = " << c.value << "°";
            break;
        case ConstraintType::Tangent:
            std::cout << "Tangent between curve " << c.entityA
                << " and curve " << c.entityB;
            break;
        case ConstraintType::Coincident:
            std::cout << "Coincident point " << c.entityA
                << " on curve " << c.entityB;
            break;
        }
        std::cout << "\n";
    }
}

// === 3D-Integration Methoden ===
bool Sketch::isClosed() const {
    if (curves_.empty()) return false;

    // Einfache Prüfung für den Anfang
    for (const auto& curve : curves_) {
        if (curve.type == CurveType::Circle || curve.type == CurveType::Rectangle) {
            return true;  // Kreise und Rechtecke sind immer geschlossen
        }
    }

    return false;
}

bool Sketch::isValidForExtrusion() const {
    // Grundlegende Prüfungen für Extrusion
    if (curves_.empty()) return false;
    if (!isClosed()) return false;

    return true;
}

bool Sketch::isValidForRevolve() const {
    // Für Rotation muss die Skizze eine geschlossene Kontur haben
    // und eine Rotationsachse definiert sein
    return isClosed() && !curves_.empty();
}

std::vector<Point> Sketch::getBoundaryPoints() const {
    std::vector<Point> boundary;

    for (const auto& curve : curves_) {
        switch (curve.type) {
        case CurveType::Line:
            boundary.push_back(curve.start);
            boundary.push_back(curve.end);
            break;

        case CurveType::Circle:
            // Diskretisiere Kreis in 36 Punkte
            for (int i = 0; i < 36; ++i) {
                double angle = i * 2.0 * M_PI / 36.0;
                boundary.push_back({
                    curve.center.x + curve.radius * std::cos(angle),
                    curve.center.y + curve.radius * std::sin(angle)
                    });
            }
            break;

        case CurveType::Arc:
        {  // Eigenen Scope für Variablendeklaration
            int segments = static_cast<int>((curve.endAngle - curve.startAngle) / 10.0);
            if (segments < 1) segments = 1;
            for (int i = 0; i <= segments; ++i) {
                double angle = curve.startAngle +
                    i * (curve.endAngle - curve.startAngle) / segments;
                double rad = angle * M_PI / 180.0;
                boundary.push_back({
                    curve.center.x + curve.radius * std::cos(rad),
                    curve.center.y + curve.radius * std::sin(rad)
                    });
            }
        }
        break;

        case CurveType::Rectangle:
            // 4 Eckpunkte im Uhrzeigersinn
            boundary.push_back(curve.start);
            boundary.push_back({ curve.corner2.x, curve.start.y });
            boundary.push_back(curve.corner2);
            boundary.push_back({ curve.start.x, curve.corner2.y });
            break;

        case CurveType::Spline:
            // Spline-Punkte hinzufügen
            for (const auto& p : curve.controlPoints) {
                boundary.push_back(p);
            }
            break;
        }
    }

    return boundary;
}

std::vector<Sketch::SimpleEdge> Sketch::getEdges() const {
    std::vector<SimpleEdge> edges;

    for (const auto& curve : curves_) {
        SimpleEdge edge;
        edge.curve = curve;

        // Bestimme Start- und Endpunkt für jede Kurve
        switch (curve.type) {
        case CurveType::Line:
            edge.start = curve.start;
            edge.end = curve.end;
            break;

        case CurveType::Circle:
            // Kreis hat keinen eindeutigen Start/Ende
            edge.start = { curve.center.x + curve.radius, curve.center.y };
            edge.end = edge.start;  // Geschlossen
            break;

        case CurveType::Arc:
        {  // Eigenen Scope öffnen
            edge.start = {
                curve.center.x + curve.radius * std::cos(curve.startAngle * M_PI / 180.0),
                curve.center.y + curve.radius * std::sin(curve.startAngle * M_PI / 180.0)
            };
            edge.end = {
                curve.center.x + curve.radius * std::cos(curve.endAngle * M_PI / 180.0),
                curve.center.y + curve.radius * std::sin(curve.endAngle * M_PI / 180.0)
            };
        }  // Scope schließen
        break;

        case CurveType::Rectangle:
            // Rechteck hat 4 Kanten, hier nur die erste
            edge.start = curve.start;
            edge.end = { curve.corner2.x, curve.start.y };
            break;

        default:
            continue;
        }

        edges.push_back(edge);
    }

    return edges;
}


double Sketch::calculateArea() const {
    double area = 0.0;

    // Einfache Berechnung für Rechtecke und Kreise
    for (const auto& curve : curves_) {
        switch (curve.type) {
        case CurveType::Rectangle:
            area += std::abs((curve.corner2.x - curve.start.x) *
                (curve.corner2.y - curve.start.y));
            break;

        case CurveType::Circle:
            area += M_PI * curve.radius * curve.radius;
            break;
        }
    }

    return area;
}

Point Sketch::calculateCentroid() const {
    if (curves_.empty()) return { 0, 0 };

    double sumX = 0, sumY = 0;
    int count = 0;

    // Vereinfachte Berechnung: Mittelwert aller Punkte
    auto points = getBoundaryPoints();
    for (const auto& p : points) {
        sumX += p.x;
        sumY += p.y;
        count++;
    }

    if (count > 0) {
        return { sumX / count, sumY / count };
    }

    return { 0, 0 };
}

void Sketch::transform(const std::function<Point(const Point&)>& transformFunc) {
    // Transformiere alle Punkte in den Kurven
    for (auto& curve : curves_) {
        switch (curve.type) {
        case CurveType::Line:
            curve.start = transformFunc(curve.start);
            curve.end = transformFunc(curve.end);
            break;

        case CurveType::Circle:
        case CurveType::Arc:
            curve.center = transformFunc(curve.center);
            break;

        case CurveType::Rectangle:
            curve.start = transformFunc(curve.start);
            curve.corner2 = transformFunc(curve.corner2);
            break;

        case CurveType::Spline:
            for (auto& p : curve.controlPoints) {
                p = transformFunc(p);
            }
            break;
        }
    }

    // Transformiere alle Punkte
    for (auto& p : points_) {
        p = transformFunc(p);
    }
}

// === Constraint Solver (Grundversion) ===
bool Sketch::solveConstraints(double tolerance, int maxIterations) {
    // Einfacher iterativer Löser für den Anfang
    for (int iter = 0; iter < maxIterations; ++iter) {
        bool changed = false;

        for (const auto& constraint : constraints_) {
            switch (constraint.type) {
            case ConstraintType::Horizontal:
                applyHorizontalConstraint(constraint.entityA);
                changed = true;
                break;

            case ConstraintType::Vertical:
                applyVerticalConstraint(constraint.entityA);
                changed = true;
                break;

            case ConstraintType::Distance:
                applyDistanceConstraint(constraint);
                changed = true;
                break;

            case ConstraintType::Angle:
                applyAngleConstraint(constraint);
                changed = true;
                break;

            default:
                break;
            }
        }

        if (!changed) {
            return true;  // Konvergiert
        }
    }

    return false;  // Maximale Iterationen erreicht
}

bool Sketch::isFullyConstrained() const {
    // Einfache Heuristik: Wenn Anzahl Constraints >= Anzahl Freiheitsgrade
    int dof = degreesOfFreedom();
    return dof <= 0;
}

int Sketch::degreesOfFreedom() const {
    // Vereinfachte Berechnung:
    // Jeder Punkt: 2 Freiheitsgrade (x, y)
    // Jede Kurve: zusätzliche Freiheitsgrade basierend auf Typ
    int dof = points_.size() * 2;

    for (const auto& curve : curves_) {
        switch (curve.type) {
        case CurveType::Line:
            dof += 4;  // Start- und Endpunkt
            break;
        case CurveType::Circle:
            dof += 3;  // Center (x,y) + Radius
            break;
        case CurveType::Arc:
            dof += 5;  // Center, Radius, Startwinkel, Endwinkel
            break;
        case CurveType::Rectangle:
            dof += 4;  // 2 Punkte
            break;
        case CurveType::Spline:
            dof += curve.controlPoints.size() * 2;
            break;
        }
    }

    // Jedes Constraint reduziert Freiheitsgrade
    dof -= constraints_.size();

    return std::max(0, dof);
}

// Hilfsmethoden für Constraint Solver
void Sketch::applyHorizontalConstraint(int curveIndex) {
    if (!validCurveIndex(curveIndex)) return;

    auto& curve = curves_[curveIndex];
    if (curve.type == CurveType::Line) {
        // Mache Linie horizontal: y-Koordinaten gleichsetzen
        curve.end.y = curve.start.y;
    }
}

void Sketch::applyVerticalConstraint(int curveIndex) {
    if (!validCurveIndex(curveIndex)) return;

    auto& curve = curves_[curveIndex];
    if (curve.type == CurveType::Line) {
        // Mache Linie vertikal: x-Koordinaten gleichsetzen
        curve.end.x = curve.start.x;
    }
}

void Sketch::applyDistanceConstraint(const Constraint& c) {
    // Einfache Implementierung: Abstand zwischen zwei Punkten
    if (c.entityA >= 0 && c.entityA < static_cast<int>(points_.size()) &&
        c.entityB >= 0 && c.entityB < static_cast<int>(points_.size())) {

        auto& p1 = points_[c.entityA];
        auto& p2 = points_[c.entityB];

        // Berechne aktuellen Vektor und skalieren auf gewünschte Länge
        double dx = p2.x - p1.x;
        double dy = p2.y - p1.y;
        double currentDist = std::sqrt(dx * dx + dy * dy);

        if (currentDist > 0) {
            double scale = c.value / currentDist;
            p2.x = p1.x + dx * scale;
            p2.y = p1.y + dy * scale;
        }
    }
}

void Sketch::applyAngleConstraint(const Constraint& c) {
    if (!validCurveIndex(c.entityA) || !validCurveIndex(c.entityB)) return;

    auto& curve1 = curves_[c.entityA];
    auto& curve2 = curves_[c.entityB];

    if (curve1.type == CurveType::Line && curve2.type == CurveType::Line) {
        // Vereinfachte Implementierung: Rotiere zweite Linie
        double currentAngle = std::atan2(curve1.end.y - curve1.start.y,
            curve1.end.x - curve1.start.x);
        double targetAngle = currentAngle + c.value * M_PI / 180.0;

        double length = std::sqrt(std::pow(curve2.end.x - curve2.start.x, 2) +
            std::pow(curve2.end.y - curve2.start.y, 2));

        curve2.end.x = curve2.start.x + length * std::cos(targetAngle);
        curve2.end.y = curve2.start.y + length * std::sin(targetAngle);
    }
}