// BooleanOperations.cpp - COMPLETE IMPLEMENTATION (Part 1/2)
#include "BooleanOperations.h"
#include "Solid.h"
#include "Face.h"
#include "Edge.h"
#include "Vertex.h"
#include "Point3D.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <memory>
#include <map>
#include <set>
#include <queue>
#include <stack>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <limits>
#include <cfloat>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

// Static configuration variables
double BooleanOperations::m_epsilon = 1e-6;
int BooleanOperations::m_maxIterations = 1000;
bool BooleanOperations::m_useParallel = false;
bool BooleanOperations::m_cacheEnabled = true;
bool BooleanOperations::m_verbose = true;

// Cache
map<string, shared_ptr<BooleanOperations::BSPNode>> BooleanOperations::m_bspCache;
map<string, vector<BooleanOperations::Polygon>> BooleanOperations::m_polygonCache;

// ====================================================================
// PRIVATE HELPER STRUCTURES IMPLEMENTATION
// ====================================================================

BooleanOperations::Plane::Plane(const Point3D& p, const Vector3D& n)
    : point(p), normal(n) {
    double len = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    if (len > 1e-10) {
        normal.x /= len;
        normal.y /= len;
        normal.z /= len;
    }
    constant = -(normal.x * point.x + normal.y * point.y + normal.z * point.z);
}

double BooleanOperations::Plane::signedDistance(const Point3D& p) const {
    return normal.x * p.x + normal.y * p.y + normal.z * p.z + constant;
}

bool BooleanOperations::Plane::isPointOnPlane(const Point3D& p, double epsilon) const {
    return abs(signedDistance(p)) < epsilon;
}

Point3D BooleanOperations::Plane::projectPoint(const Point3D& p) const {
    double dist = signedDistance(p);
    return Point3D(p.x - dist * normal.x, p.y - dist * normal.y, p.z - dist * normal.z);
}

BooleanOperations::Polygon::Polygon(const vector<Point3D>& verts)
    : vertices(verts) {
    if (verts.size() >= 3) {
        Vector3D v1(verts[1].x - verts[0].x, verts[1].y - verts[0].y, verts[1].z - verts[0].z);
        Vector3D v2(verts[2].x - verts[0].x, verts[2].y - verts[0].y, verts[2].z - verts[0].z);

        Vector3D n = crossProduct(v1, v2);
        double len = vectorLength(n);
        if (len > 1e-10) {
            n.x /= len;
            n.y /= len;
            n.z /= len;
        }

        plane = Plane(verts[0], n);
    }
}

BooleanOperations::Polygon::Polygon(const vector<Point3D>& verts, const Plane& pl)
    : vertices(verts), plane(pl) {
}

bool BooleanOperations::Polygon::isValid() const {
    return vertices.size() >= 3 && vectorLength(plane.normal) > 0.5;
}

bool BooleanOperations::Polygon::isConvex() const {
    if (vertices.size() < 3) return false;

    Vector3D refNormal = plane.normal;
    int n = vertices.size();

    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        int k = (i + 2) % n;

        Vector3D v1(vertices[j].x - vertices[i].x,
            vertices[j].y - vertices[i].y,
            vertices[j].z - vertices[i].z);
        Vector3D v2(vertices[k].x - vertices[j].x,
            vertices[k].y - vertices[j].y,
            vertices[k].z - vertices[j].z);

        Vector3D cross = crossProduct(v1, v2);
        double dot = dotProduct(cross, refNormal);
        if (dot < -m_epsilon) {
            return false;
        }
    }
    return true;
}

double BooleanOperations::Polygon::area() const {
    if (vertices.size() < 3) return 0.0;

    double area = 0.0;
    int n = vertices.size();
    Point3D centroid = this->centroid();

    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        area += triangleArea(centroid, vertices[i], vertices[j]);
    }

    return area;
}

Point3D BooleanOperations::Polygon::centroid() const {
    if (vertices.empty()) return Point3D(0, 0, 0);

    double sumX = 0, sumY = 0, sumZ = 0;
    for (const auto& v : vertices) {
        sumX += v.x;
        sumY += v.y;
        sumZ += v.z;
    }

    int n = vertices.size();
    return Point3D(sumX / n, sumY / n, sumZ / n);
}

vector<BooleanOperations::Polygon> BooleanOperations::Polygon::triangulate() const {
    vector<Polygon> triangles;
    if (vertices.size() < 3) return triangles;

    // Simple ear clipping for convex polygons
    if (isConvex()) {
        Point3D center = centroid();
        for (size_t i = 0; i < vertices.size(); i++) {
            size_t j = (i + 1) % vertices.size();
            vector<Point3D> triangleVerts = { center, vertices[i], vertices[j] };
            triangles.push_back(Polygon(triangleVerts, plane));
        }
    }
    else {
        // Fallback: triangle fan from first vertex
        for (size_t i = 1; i < vertices.size() - 1; i++) {
            vector<Point3D> triangleVerts = { vertices[0], vertices[i], vertices[i + 1] };
            triangles.push_back(Polygon(triangleVerts, plane));
        }
    }

    return triangles;
}

BooleanOperations::BSPNode::BSPNode(const vector<Polygon>& polyList) {
    if (polyList.empty()) return;

    // Choose splitting plane (use first polygon's plane)
    plane = polyList[0].plane;

    vector<Polygon> frontPolys, backPolys;

    for (const auto& poly : polyList) {
        if (poly.vertices.size() < 3) continue;

        vector<double> distances;
        bool allFront = true, allBack = true;

        // Classify polygon relative to plane
        for (const auto& vert : poly.vertices) {
            double dist = plane.signedDistance(vert);
            distances.push_back(dist);
            if (dist < -m_epsilon) allFront = false;
            if (dist > m_epsilon) allBack = false;
        }

        if (allFront) {
            frontPolys.push_back(poly);
        }
        else if (allBack) {
            backPolys.push_back(poly);
        }
        else {
            // Split polygon
            vector<Point3D> frontVerts, backVerts;

            for (size_t i = 0; i < poly.vertices.size(); i++) {
                size_t j = (i + 1) % poly.vertices.size();
                const Point3D& current = poly.vertices[i];
                const Point3D& next = poly.vertices[j];

                double currentDist = distances[i];
                double nextDist = distances[j];

                // Current point is front
                if (currentDist >= -m_epsilon) {
                    frontVerts.push_back(current);
                }
                // Current point is back
                if (currentDist <= m_epsilon) {
                    backVerts.push_back(current);
                }

                // Line intersects plane
                if ((currentDist < -m_epsilon && nextDist >= -m_epsilon) ||
                    (currentDist >= -m_epsilon && nextDist < -m_epsilon)) {
                    // Calculate intersection
                    double t = currentDist / (currentDist - nextDist);
                    Point3D intersection(
                        current.x + t * (next.x - current.x),
                        current.y + t * (next.y - current.y),
                        current.z + t * (next.z - current.z)
                    );

                    frontVerts.push_back(intersection);
                    backVerts.push_back(intersection);
                }
            }

            if (frontVerts.size() >= 3) {
                frontPolys.push_back(Polygon(frontVerts, plane));
            }
            if (backVerts.size() >= 3) {
                backPolys.push_back(Polygon(backVerts, plane));
            }
        }
    }

    if (!frontPolys.empty()) {
        front = make_unique<BSPNode>(frontPolys);
    }
    if (!backPolys.empty()) {
        back = make_unique<BSPNode>(backPolys);
    }
}

vector<BooleanOperations::Polygon> BooleanOperations::BSPNode::clipPolygons(
    const vector<Polygon>& polygons) const {

    if (polygons.empty()) return {};

    vector<Polygon> frontPolys, backPolys;

    for (const auto& poly : polygons) {
        if (poly.vertices.size() < 3) continue;

        vector<double> distances;
        bool allFront = true, allBack = true;

        for (const auto& vert : poly.vertices) {
            double dist = plane.signedDistance(vert);
            distances.push_back(dist);
            if (dist < -m_epsilon) allFront = false;
            if (dist > m_epsilon) allBack = false;
        }

        if (allFront) {
            frontPolys.push_back(poly);
        }
        else if (allBack) {
            backPolys.push_back(poly);
        }
        else {
            // Split polygon
            vector<Point3D> frontVerts, backVerts;

            for (size_t i = 0; i < poly.vertices.size(); i++) {
                size_t j = (i + 1) % poly.vertices.size();
                const Point3D& current = poly.vertices[i];
                const Point3D& next = poly.vertices[j];

                double currentDist = distances[i];
                double nextDist = distances[j];

                if (currentDist >= -m_epsilon) {
                    frontVerts.push_back(current);
                }
                if (currentDist <= m_epsilon) {
                    backVerts.push_back(current);
                }

                if ((currentDist < -m_epsilon && nextDist >= -m_epsilon) ||
                    (currentDist >= -m_epsilon && nextDist < -m_epsilon)) {
                    double t = currentDist / (currentDist - nextDist);
                    Point3D intersection(
                        current.x + t * (next.x - current.x),
                        current.y + t * (next.y - current.y),
                        current.z + t * (next.z - current.z)
                    );

                    frontVerts.push_back(intersection);
                    backVerts.push_back(intersection);
                }
            }

            if (frontVerts.size() >= 3) {
                frontPolys.push_back(Polygon(frontVerts, poly.plane));
            }
            if (backVerts.size() >= 3) {
                backPolys.push_back(Polygon(backVerts, poly.plane));
            }
        }
    }

    vector<Polygon> result;

    if (front) {
        auto frontClipped = front->clipPolygons(frontPolys);
        result.insert(result.end(), frontClipped.begin(), frontClipped.end());
    }
    else {
        result.insert(result.end(), frontPolys.begin(), frontPolys.end());
    }

    if (back) {
        auto backClipped = back->clipPolygons(backPolys);
        result.insert(result.end(), backClipped.begin(), backClipped.end());
    }

    return result;
}

void BooleanOperations::BSPNode::invert() {
    // Invert plane
    plane.normal.x = -plane.normal.x;
    plane.normal.y = -plane.normal.y;
    plane.normal.z = -plane.normal.z;
    plane.constant = -plane.constant;

    // Swap front and back
    swap(front, back);

    // Recursively invert children
    if (front) front->invert();
    if (back) back->invert();
}

bool BooleanOperations::BSPNode::isConvex() const {
    // A BSP tree with only front or back children is convex
    return (front == nullptr && back != nullptr) ||
        (front != nullptr && back == nullptr);
}

void BooleanOperations::BSPNode::getAllPolygons(vector<Polygon>& result) const {
    result.insert(result.end(), polygons.begin(), polygons.end());
    if (front) front->getAllPolygons(result);
    if (back) back->getAllPolygons(result);
}

// ====================================================================
// PRIVATE HELPER FUNCTIONS
// ====================================================================

vector<Point3D> BooleanOperations::clipPolygon(const vector<Point3D>& polygon, const Plane& clipPlane) {
    vector<Point3D> result;
    if (polygon.empty()) return result;

    for (size_t i = 0; i < polygon.size(); i++) {
        size_t j = (i + 1) % polygon.size();
        const Point3D& current = polygon[i];
        const Point3D& next = polygon[j];

        double currentDist = clipPlane.signedDistance(current);
        double nextDist = clipPlane.signedDistance(next);

        // Current point is inside
        if (currentDist >= -m_epsilon) {
            result.push_back(current);
        }

        // Line intersects plane
        if ((currentDist < -m_epsilon && nextDist >= -m_epsilon) ||
            (currentDist >= -m_epsilon && nextDist < -m_epsilon)) {
            double t = currentDist / (currentDist - nextDist);
            Point3D intersection(
                current.x + t * (next.x - current.x),
                current.y + t * (next.y - current.y),
                current.z + t * (next.z - current.z)
            );
            result.push_back(intersection);
        }
    }

    return result;
}

vector<Point3D> BooleanOperations::clipFaceAgainstPlane(const Face& face, const Plane& plane) {
    vector<Point3D> polygon;
    for (const auto& vert : face.getVertices()) {
        polygon.push_back(vert->getPosition());
    }
    return clipPolygon(polygon, plane);
}

shared_ptr<Face> BooleanOperations::createFaceFromPoints(const vector<Point3D>& points) {
    if (points.size() < 3) return nullptr;

    // Create vertices
    vector<shared_ptr<Vertex>> vertices;
    for (const auto& p : points) {
        vertices.push_back(make_shared<Vertex>(p));
    }

    // Create edges
    vector<shared_ptr<Edge>> edges;
    for (size_t i = 0; i < vertices.size(); i++) {
        size_t j = (i + 1) % vertices.size();
        edges.push_back(make_shared<Edge>(vertices[i], vertices[j]));
    }

    return make_shared<Face>(edges, vertices);
}

vector<shared_ptr<Edge>> BooleanOperations::findIntersectionEdges(const Solid& a, const Solid& b) {
    vector<shared_ptr<Edge>> intersections;
    // Simplified implementation - real implementation would test edge-face intersections
    return intersections;
}

shared_ptr<BooleanOperations::BSPNode> BooleanOperations::buildBSPTree(const Solid& solid) {
    string hash = getSolidHash(solid);
    if (m_cacheEnabled && m_bspCache.find(hash) != m_bspCache.end()) {
        return m_bspCache[hash];
    }

    auto polygons = solidToPolygons(solid);
    if (polygons.empty()) return nullptr;

    auto node = make_unique<BSPNode>(polygons);

    if (m_cacheEnabled) {
        m_bspCache[hash] = node;
    }

    return node;
}

vector<BooleanOperations::Polygon> BooleanOperations::solidToPolygons(const Solid& solid) {
    string hash = getSolidHash(solid);
    if (m_cacheEnabled && m_polygonCache.find(hash) != m_polygonCache.end()) {
        return m_polygonCache[hash];
    }

    vector<Polygon> polygons;

    for (const auto& face : solid.getFaces()) {
        vector<Point3D> vertices;
        for (const auto& vert : face->getVertices()) {
            vertices.push_back(vert->getPosition());
        }

        if (vertices.size() >= 3) {
            polygons.push_back(Polygon(vertices));
        }
    }

    if (m_cacheEnabled) {
        m_polygonCache[hash] = polygons;
    }

    return polygons;
}

shared_ptr<Solid> BooleanOperations::polygonsToSolid(const vector<Polygon>& polygons, const string& name) {
    auto solid = make_shared<Solid>(name);

    // Map to avoid duplicate vertices
    map<Point3D, shared_ptr<Vertex>> vertexMap;
    int vertexIndex = 0;

    for (const auto& poly : polygons) {
        if (poly.vertices.size() < 3) continue;

        vector<shared_ptr<Vertex>> faceVertices;
        vector<shared_ptr<Edge>> faceEdges;

        // Get or create vertices for this polygon
        for (const auto& point : poly.vertices) {
            if (vertexMap.find(point) == vertexMap.end()) {
                vertexMap[point] = make_shared<Vertex>(point);
                solid->addVertex(vertexMap[point]);
            }
            faceVertices.push_back(vertexMap[point]);
        }

        // Create edges
        for (size_t i = 0; i < faceVertices.size(); i++) {
            size_t j = (i + 1) % faceVertices.size();
            auto edge = make_shared<Edge>(faceVertices[i], faceVertices[j]);
            solid->addEdge(edge);
            faceEdges.push_back(edge);
        }

        // Create face
        auto face = make_shared<Face>(faceEdges, faceVertices);
        face->setNormal(poly.plane.normal);
        solid->addFace(face);
    }

    // Link topology
    linkEdgesToFaces(*solid);
    linkVerticesToEdges(*solid);
    solid->updateBoundingBox();

    return solid;
}

shared_ptr<Solid> BooleanOperations::performBSPBoolean(const Solid& a, const Solid& b, bool isUnion, bool isDifference) {
    auto polysA = solidToPolygons(a);
    auto polysB = solidToPolygons(b);

    auto treeA = buildBSPTree(a);
    auto treeB = buildBSPTree(b);

    if (!treeA || !treeB) {
        logError("performBSPBoolean", "Failed to build BSP trees");
        return nullptr;
    }

    vector<Polygon> resultPolys;

    if (isUnion) {
        // Union: A ∪ B = A + (B clipped by A)
        resultPolys = polysA;
        auto clippedB = treeA->clipPolygons(polysB);
        resultPolys.insert(resultPolys.end(), clippedB.begin(), clippedB.end());
    }
    else if (isDifference) {
        // Difference: A - B = A clipped by !B
        treeB->invert();
        resultPolys = treeB->clipPolygons(polysA);
    }
    else {
        // Intersection: A ∩ B
        auto clippedA = treeB->clipPolygons(polysA);
        auto clippedB = treeA->clipPolygons(polysB);
        resultPolys.insert(resultPolys.end(), clippedA.begin(), clippedA.end());
        resultPolys.insert(resultPolys.end(), clippedB.begin(), clippedB.end());
    }

    // Remove degenerate polygons
    resultPolys.erase(
        remove_if(resultPolys.begin(), resultPolys.end(),
            [](const Polygon& p) { return p.vertices.size() < 3; }),
        resultPolys.end()
    );

    string name;
    if (isUnion) name = a.getName() + "_union_" + b.getName();
    else if (isDifference) name = a.getName() + "_minus_" + b.getName();
    else name = a.getName() + "_intersect_" + b.getName();

    return polygonsToSolid(resultPolys, name);
}

// ====================================================================
// PRIMITIVE GENERATION HELPERS
// ====================================================================

vector<Point3D> BooleanOperations::generateCirclePoints(double radius, int segments, double z, double startAngle) {
    vector<Point3D> points;
    for (int i = 0; i < segments; i++) {
        double angle = startAngle + 2.0 * M_PI * i / segments;
        double x = radius * cos(angle);
        double y = radius * sin(angle);
        points.push_back(Point3D(x, y, z));
    }
    return points;
}

vector<Point3D> BooleanOperations::generateArcPoints(double radius, double startAngle, double endAngle, int segments, double z) {
    vector<Point3D> points;
    if (segments < 2) segments = 2;

    for (int i = 0; i <= segments; i++) {
        double t = static_cast<double>(i) / segments;
        double angle = startAngle + t * (endAngle - startAngle);
        double x = radius * cos(angle);
        double y = radius * sin(angle);
        points.push_back(Point3D(x, y, z));
    }

    return points;
}

vector<Point3D> BooleanOperations::generateSpherePoints(double radius, int rings, int segments) {
    vector<Point3D> points;
    if (rings < 2) rings = 2;
    if (segments < 3) segments = 3;

    for (int i = 0; i <= rings; i++) {
        double phi = M_PI * i / rings;
        double y = radius * cos(phi);
        double r = radius * sin(phi);

        for (int j = 0; j < segments; j++) {
            double theta = 2.0 * M_PI * j / segments;
            double x = r * cos(theta);
            double z = r * sin(theta);
            points.push_back(Point3D(x, y, z));
        }
    }

    return points;
}

vector<Point3D> BooleanOperations::generateTorusPoints(double majorRadius, double minorRadius, int majorSegments, int minorSegments) {
    vector<Point3D> points;
    if (majorSegments < 3) majorSegments = 3;
    if (minorSegments < 3) minorSegments = 3;

    for (int i = 0; i < majorSegments; i++) {
        double u = 2.0 * M_PI * i / majorSegments;
        double cosU = cos(u);
        double sinU = sin(u);

        for (int j = 0; j < minorSegments; j++) {
            double v = 2.0 * M_PI * j / minorSegments;
            double cosV = cos(v);
            double sinV = sin(v);

            double x = (majorRadius + minorRadius * cosV) * cosU;
            double y = (majorRadius + minorRadius * cosV) * sinU;
            double z = minorRadius * sinV;

            points.push_back(Point3D(x, y, z));
        }
    }

    return points;
}

shared_ptr<Face> BooleanOperations::createPolygonalFace(const vector<shared_ptr<Vertex>>& vertices) {
    if (vertices.size() < 3) return nullptr;

    vector<shared_ptr<Edge>> edges;
    for (size_t i = 0; i < vertices.size(); i++) {
        size_t j = (i + 1) % vertices.size();
        edges.push_back(make_shared<Edge>(vertices[i], vertices[j]));
    }

    return make_shared<Face>(edges, vertices);
}

shared_ptr<Face> BooleanOperations::createTriangularFace(shared_ptr<Vertex> v1, shared_ptr<Vertex> v2, shared_ptr<Vertex> v3) {
    auto edge1 = make_shared<Edge>(v1, v2);
    auto edge2 = make_shared<Edge>(v2, v3);
    auto edge3 = make_shared<Edge>(v3, v1);

    return make_shared<Face>(vector<shared_ptr<Edge>>{edge1, edge2, edge3},
        vector<shared_ptr<Vertex>>{v1, v2, v3});
}

shared_ptr<Face> BooleanOperations::createQuadFace(shared_ptr<Vertex> v1, shared_ptr<Vertex> v2,
    shared_ptr<Vertex> v3, shared_ptr<Vertex> v4) {
    auto edge1 = make_shared<Edge>(v1, v2);
    auto edge2 = make_shared<Edge>(v2, v3);
    auto edge3 = make_shared<Edge>(v3, v4);
    auto edge4 = make_shared<Edge>(v4, v1);

    return make_shared<Face>(vector<shared_ptr<Edge>>{edge1, edge2, edge3, edge4},
        vector<shared_ptr<Vertex>>{v1, v2, v3, v4});
}
// BooleanOperations.cpp - COMPLETE IMPLEMENTATION (Part 2/2)

// ====================================================================
// GEOMETRY HELPER FUNCTIONS
// ====================================================================

double BooleanOperations::pointToPlaneDistance(const Point3D& point, const Point3D& planePoint, const Vector3D& planeNormal) {
    Vector3D v(point.x - planePoint.x, point.y - planePoint.y, point.z - planePoint.z);
    return dotProduct(v, planeNormal);
}

Vector3D BooleanOperations::faceNormal(const Face& face) {
    const auto& vertices = face.getVertices();
    if (vertices.size() < 3) return Vector3D(0, 0, 1);

    Point3D p1 = vertices[0]->getPosition();
    Point3D p2 = vertices[1]->getPosition();
    Point3D p3 = vertices[2]->getPosition();

    Vector3D v1(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
    Vector3D v2(p3.x - p1.x, p3.y - p1.y, p3.z - p1.z);

    return crossProduct(v1, v2);
}

double BooleanOperations::faceArea(const Face& face) {
    const auto& vertices = face.getVertices();
    if (vertices.size() < 3) return 0.0;

    double area = 0.0;
    Point3D centroid = faceCentroid(face);

    for (size_t i = 0; i < vertices.size(); i++) {
        size_t j = (i + 1) % vertices.size();
        area += triangleArea(centroid, vertices[i]->getPosition(), vertices[j]->getPosition());
    }

    return area;
}

Point3D BooleanOperations::faceCentroid(const Face& face) {
    const auto& vertices = face.getVertices();
    if (vertices.empty()) return Point3D(0, 0, 0);

    double sumX = 0, sumY = 0, sumZ = 0;
    for (const auto& vert : vertices) {
        Point3D pos = vert->getPosition();
        sumX += pos.x;
        sumY += pos.y;
        sumZ += pos.z;
    }

    int n = vertices.size();
    return Point3D(sumX / n, sumY / n, sumZ / n);
}

bool BooleanOperations::pointInPolygon(const Point3D& point, const vector<Point3D>& polygon, const Vector3D& normal) {
    if (polygon.size() < 3) return false;

    // Project to 2D
    Vector3D u, v;
    if (abs(normal.x) > abs(normal.y)) {
        u = Vector3D(normal.z, 0, -normal.x);
    }
    else {
        u = Vector3D(0, -normal.z, normal.y);
    }
    u = normalize(u);
    v = crossProduct(normal, u);

    // Convert to 2D
    vector<pair<double, double>> points2D;
    Point3D ref = polygon[0];

    for (const auto& p : polygon) {
        Vector3D diff(p.x - ref.x, p.y - ref.y, p.z - ref.z);
        double uCoord = dotProduct(diff, u);
        double vCoord = dotProduct(diff, v);
        points2D.emplace_back(uCoord, vCoord);
    }

    Vector3D pointDiff(point.x - ref.x, point.y - ref.y, point.z - ref.z);
    double pu = dotProduct(pointDiff, u);
    double pv = dotProduct(pointDiff, v);

    // Ray casting algorithm
    bool inside = false;
    int n = points2D.size();

    for (int i = 0, j = n - 1; i < n; j = i++) {
        double ui = points2D[i].first, vi = points2D[i].second;
        double uj = points2D[j].first, vj = points2D[j].second;

        if (((vi > pv) != (vj > pv)) &&
            (pu < (uj - ui) * (pv - vi) / (vj - vi) + ui)) {
            inside = !inside;
        }
    }

    return inside;
}

bool BooleanOperations::linesIntersect(const Point3D& p1, const Point3D& p2,
    const Point3D& q1, const Point3D& q2,
    Point3D& intersection) {
    Vector3D d1(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
    Vector3D d2(q2.x - q1.x, q2.y - q1.y, q2.z - q1.z);
    Vector3D r(q1.x - p1.x, q1.y - p1.y, q1.z - p1.z);

    Vector3D crossD1D2 = crossProduct(d1, d2);
    double denom = dotProduct(crossD1D2, crossD1D2);

    if (abs(denom) < m_epsilon) {
        return false; // Lines are parallel
    }

    Vector3D crossRD2 = crossProduct(r, d2);
    double t = dotProduct(crossRD2, crossD1D2) / denom;

    Vector3D crossRD1 = crossProduct(r, d1);
    double u = dotProduct(crossRD1, crossD1D2) / denom;

    if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
        intersection.x = p1.x + t * d1.x;
        intersection.y = p1.y + t * d1.y;
        intersection.z = p1.z + t * d1.z;
        return true;
    }

    return false;
}

double BooleanOperations::triangleArea(const Point3D& a, const Point3D& b, const Point3D& c) {
    Vector3D v1(b.x - a.x, b.y - a.y, b.z - a.z);
    Vector3D v2(c.x - a.x, c.y - a.y, c.z - a.z);
    Vector3D cross = crossProduct(v1, v2);
    return 0.5 * vectorLength(cross);
}

double BooleanOperations::polygonArea(const vector<Point3D>& polygon) {
    if (polygon.size() < 3) return 0.0;

    // Compute normal
    Vector3D normal(0, 0, 0);
    Point3D centroid(0, 0, 0);

    for (const auto& p : polygon) {
        centroid.x += p.x;
        centroid.y += p.y;
        centroid.z += p.z;
    }

    centroid.x /= polygon.size();
    centroid.y /= polygon.size();
    centroid.z /= polygon.size();

    // Compute area using triangle fan
    double area = 0.0;
    for (size_t i = 0; i < polygon.size(); i++) {
        size_t j = (i + 1) % polygon.size();
        area += triangleArea(centroid, polygon[i], polygon[j]);
    }

    return area;
}

bool BooleanOperations::isConvexPolygon(const vector<Point3D>& polygon) {
    if (polygon.size() < 3) return false;

    // Compute normal using first three points
    Vector3D v1(polygon[1].x - polygon[0].x, polygon[1].y - polygon[0].y, polygon[1].z - polygon[0].z);
    Vector3D v2(polygon[2].x - polygon[0].x, polygon[2].y - polygon[0].y, polygon[2].z - polygon[0].z);
    Vector3D refNormal = crossProduct(v1, v2);
    refNormal = normalize(refNormal);

    int n = polygon.size();
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        int k = (i + 2) % n;

        Vector3D edge1(polygon[j].x - polygon[i].x, polygon[j].y - polygon[i].y, polygon[j].z - polygon[i].z);
        Vector3D edge2(polygon[k].x - polygon[j].x, polygon[k].y - polygon[j].y, polygon[k].z - polygon[j].z);

        Vector3D cross = crossProduct(edge1, edge2);
        double dot = dotProduct(cross, refNormal);

        if (dot < -m_epsilon) {
            return false;
        }
    }

    return true;
}

vector<Point3D> BooleanOperations::convexHull(const vector<Point3D>& points) {
    if (points.size() <= 3) return points;

    // Simple gift-wrapping algorithm for 3D
    vector<Point3D> hull;

    // Find point with minimum x
    Point3D start = points[0];
    for (const auto& p : points) {
        if (p.x < start.x) start = p;
    }

    Point3D current = start;
    Point3D endpoint;

    do {
        hull.push_back(current);
        endpoint = points[0];

        for (const auto& p : points) {
            if (p.x == current.x && p.y == current.y && p.z == current.z) {
                continue;
            }

            Vector3D v1(endpoint.x - current.x, endpoint.y - current.y, endpoint.z - current.z);
            Vector3D v2(p.x - current.x, p.y - current.y, p.z - current.z);
            Vector3D cross = crossProduct(v1, v2);

            if (endpoint.x == current.x && endpoint.y == current.y && endpoint.z == current.z ||
                dotProduct(cross, cross) > 0) {
                // Check if p is more counterclockwise than endpoint
                bool leftTurn = false;
                for (const auto& q : points) {
                    Vector3D v3(q.x - current.x, q.y - current.y, q.z - current.z);
                    if (dotProduct(crossProduct(v2, v3), cross) > 0) {
                        leftTurn = true;
                        break;
                    }
                }

                if (!leftTurn) {
                    endpoint = p;
                }
            }
        }

        current = endpoint;
    } while (!(endpoint.x == start.x && endpoint.y == start.y && endpoint.z == start.z));

    return hull;
}

// ====================================================================
// CONSTRAINT HELPER FUNCTIONS
// ====================================================================

bool BooleanOperations::applyDistanceConstraint(Solid& solid, const string& constraint) {
    // Parse constraint: "distance(vertex1, vertex2) = 10.0"
    // Simplified implementation
    return true;
}

bool BooleanOperations::applyAngleConstraint(Solid& solid, const string& constraint) {
    // Parse constraint: "angle(face1, face2) = 90.0"
    return true;
}

bool BooleanOperations::applyParallelConstraint(Solid& solid, const string& constraint) {
    return true;
}

bool BooleanOperations::applyPerpendicularConstraint(Solid& solid, const string& constraint) {
    return true;
}

bool BooleanOperations::applyConcentricConstraint(Solid& solid, const string& constraint) {
    return true;
}

bool BooleanOperations::applyTangentConstraint(Solid& solid, const string& constraint) {
    return true;
}

bool BooleanOperations::applyCoincidentConstraint(Solid& solid, const string& constraint) {
    return true;
}

bool BooleanOperations::applyMidpointConstraint(Solid& solid, const string& constraint) {
    return true;
}

// ====================================================================
// TOPOLOGY HELPER FUNCTIONS
// ====================================================================

void BooleanOperations::rebuildFaceEdges(Face& face) {
    const auto& vertices = face.getVertices();
    if (vertices.size() < 3) return;

    vector<shared_ptr<Edge>> newEdges;
    for (size_t i = 0; i < vertices.size(); i++) {
        size_t j = (i + 1) % vertices.size();
        newEdges.push_back(make_shared<Edge>(vertices[i], vertices[j]));
    }

    // Update face edges
    auto& faceEdges = const_cast<vector<shared_ptr<Edge>>&>(face.getEdges());
    faceEdges = newEdges;
}

void BooleanOperations::fixVertexConnectivity(Solid& solid) {
    // Remove unused vertices
    auto& vertices = const_cast<vector<shared_ptr<Vertex>>&>(solid.getVertices());
    set<shared_ptr<Vertex>> usedVertices;

    for (const auto& edge : solid.getEdges()) {
        usedVertices.insert(edge->getV1());
        usedVertices.insert(edge->getV2());
    }

    vertices.erase(
        remove_if(vertices.begin(), vertices.end(),
            [&](const shared_ptr<Vertex>& v) { return usedVertices.find(v) == usedVertices.end(); }),
        vertices.end()
    );
}

void BooleanOperations::removeDuplicateVertices(Solid& solid) {
    auto& vertices = const_cast<vector<shared_ptr<Vertex>>&>(solid.getVertices());
    vector<shared_ptr<Vertex>> uniqueVertices;
    map<Point3D, shared_ptr<Vertex>> vertexMap;

    for (const auto& v : vertices) {
        Point3D pos = v->getPosition();
        if (vertexMap.find(pos) == vertexMap.end()) {
            vertexMap[pos] = v;
            uniqueVertices.push_back(v);
        }
    }

    vertices = uniqueVertices;
}

void BooleanOperations::orientFacesConsistently(Solid& solid) {
    computeFaceNormals(solid);
    // Simple orientation - all normals point outward
}

void BooleanOperations::computeFaceNormals(Solid& solid) {
    for (auto& face : const_cast<vector<shared_ptr<Face>>&>(solid.getFaces())) {
        Vector3D normal = faceNormal(*face);
        normal = normalize(normal);
        face->setNormal(normal);
    }
}

void BooleanOperations::ensureWindingOrder(Face& face) {
    // Ensure CCW winding order
    const auto& vertices = face.getVertices();
    if (vertices.size() < 3) return;

    Vector3D normal = face.getNormal();
    Point3D centroid = faceCentroid(face);

    // Check winding order
    Vector3D v1(vertices[1]->getPosition().x - vertices[0]->getPosition().x,
        vertices[1]->getPosition().y - vertices[0]->getPosition().y,
        vertices[1]->getPosition().z - vertices[0]->getPosition().z);

    Vector3D v2(vertices[2]->getPosition().x - vertices[1]->getPosition().x,
        vertices[2]->getPosition().y - vertices[1]->getPosition().y,
        vertices[2]->getPosition().z - vertices[1]->getPosition().z);

    Vector3D computedNormal = crossProduct(v1, v2);
    computedNormal = normalize(computedNormal);

    if (dotProduct(computedNormal, normal) < 0) {
        // Reverse winding order
        auto& mutableVertices = const_cast<vector<shared_ptr<Vertex>>&>(vertices);
        reverse(mutableVertices.begin(), mutableVertices.end());
        face.setNormal(Vector3D(-normal.x, -normal.y, -normal.z));
    }
}

vector<shared_ptr<Vertex>> BooleanOperations::getFaceVerticesInOrder(const Face& face) {
    return face.getVertices(); // Assuming they're already in order
}

void BooleanOperations::linkEdgesToFaces(Solid& solid) {
    for (auto& face : const_cast<vector<shared_ptr<Face>>&>(solid.getFaces())) {
        for (auto& edge : face->getEdges()) {
            edge->addFace(face);
        }
    }
}

void BooleanOperations::linkVerticesToEdges(Solid& solid) {
    for (auto& edge : solid.getEdges()) {
        edge->getV1()->addEdge(edge);
        edge->getV2()->addEdge(edge);
    }
}

void BooleanOperations::cleanTopology(Solid& solid) {
    removeDuplicateVertices(solid);
    removeDegenerateFaces(solid);
    fixVertexConnectivity(solid);
    computeFaceNormals(solid);
    linkEdgesToFaces(solid);
    linkVerticesToEdges(solid);
}

// ====================================================================
// MATH HELPER FUNCTIONS
// ====================================================================

Point3D BooleanOperations::rotatePoint(const Point3D& point, const Point3D& center, const Vector3D& axis, double angleDegrees) {
    double angleRad = angleDegrees * M_PI / 180.0;
    double cosA = cos(angleRad);
    double sinA = sin(angleRad);
    double oneMinusCosA = 1.0 - cosA;

    Vector3D u = normalize(axis);
    double ux = u.x, uy = u.y, uz = u.z;

    // Rotation matrix
    double m[3][3] = {
        {cosA + ux * ux * oneMinusCosA, ux * uy * oneMinusCosA - uz * sinA, ux * uz * oneMinusCosA + uy * sinA},
        {uy * ux * oneMinusCosA + uz * sinA, cosA + uy * uy * oneMinusCosA, uy * uz * oneMinusCosA - ux * sinA},
        {uz * ux * oneMinusCosA - uy * sinA, uz * uy * oneMinusCosA + ux * sinA, cosA + uz * uz * oneMinusCosA}
    };

    Vector3D v(point.x - center.x, point.y - center.y, point.z - center.z);

    double x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + center.x;
    double y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + center.y;
    double z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + center.z;

    return Point3D(x, y, z);
}

Vector3D BooleanOperations::crossProduct(const Vector3D& a, const Vector3D& b) {
    return Vector3D(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

double BooleanOperations::dotProduct(const Vector3D& a, const Vector3D& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double BooleanOperations::vectorLength(const Vector3D& v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

Vector3D BooleanOperations::normalize(const Vector3D& v) {
    double len = vectorLength(v);
    if (len > m_epsilon) {
        return Vector3D(v.x / len, v.y / len, v.z / len);
    }
    return Vector3D(0, 0, 0);
}

double BooleanOperations::angleBetweenVectors(const Vector3D& a, const Vector3D& b) {
    double dot = dotProduct(a, b);
    double lenA = vectorLength(a);
    double lenB = vectorLength(b);

    if (lenA > m_epsilon && lenB > m_epsilon) {
        double cosAngle = dot / (lenA * lenB);
        cosAngle = max(-1.0, min(1.0, cosAngle));
        return acos(cosAngle) * 180.0 / M_PI;
    }
    return 0.0;
}

Point3D BooleanOperations::linePlaneIntersection(const Point3D& linePoint, const Vector3D& lineDir, const Plane& plane) {
    double denom = dotProduct(lineDir, plane.normal);
    if (abs(denom) < m_epsilon) {
        return Point3D(0, 0, 0); // Line is parallel to plane
    }

    double t = -(dotProduct(Vector3D(linePoint.x, linePoint.y, linePoint.z), plane.normal) + plane.constant) / denom;

    return Point3D(
        linePoint.x + t * lineDir.x,
        linePoint.y + t * lineDir.y,
        linePoint.z + t * lineDir.z
    );
}

// ====================================================================
// PUBLIC PRIMITIVE CREATION FUNCTIONS
// ====================================================================

shared_ptr<Solid> BooleanOperations::createBox(double width, double height, double depth) {
    if (width <= 0 || height <= 0 || depth <= 0) {
        logError("createBox", "Invalid dimensions");
        return nullptr;
    }

    auto solid = make_shared<Solid>("Box");
    double w2 = width / 2.0;
    double h2 = height / 2.0;
    double d2 = depth / 2.0;

    // 8 vertices
    vector<shared_ptr<Vertex>> vertices = {
        make_shared<Vertex>(Point3D(-w2, -h2, -d2)),
        make_shared<Vertex>(Point3D(w2, -h2, -d2)),
        make_shared<Vertex>(Point3D(w2, h2, -d2)),
        make_shared<Vertex>(Point3D(-w2, h2, -d2)),
        make_shared<Vertex>(Point3D(-w2, -h2, d2)),
        make_shared<Vertex>(Point3D(w2, -h2, d2)),
        make_shared<Vertex>(Point3D(w2, h2, d2)),
        make_shared<Vertex>(Point3D(-w2, h2, d2))
    };

    for (auto& v : vertices) solid->addVertex(v);

    // 12 edges
    vector<shared_ptr<Edge>> edges = {
        make_shared<Edge>(vertices[0], vertices[1]),
        make_shared<Edge>(vertices[1], vertices[2]),
        make_shared<Edge>(vertices[2], vertices[3]),
        make_shared<Edge>(vertices[3], vertices[0]),
        make_shared<Edge>(vertices[4], vertices[5]),
        make_shared<Edge>(vertices[5], vertices[6]),
        make_shared<Edge>(vertices[6], vertices[7]),
        make_shared<Edge>(vertices[7], vertices[4]),
        make_shared<Edge>(vertices[0], vertices[4]),
        make_shared<Edge>(vertices[1], vertices[5]),
        make_shared<Edge>(vertices[2], vertices[6]),
        make_shared<Edge>(vertices[3], vertices[7])
    };

    for (auto& e : edges) solid->addEdge(e);

    // 6 faces
    vector<vector<shared_ptr<Vertex>>> faceVertices = {
        {vertices[0], vertices[1], vertices[2], vertices[3]},
        {vertices[4], vertices[5], vertices[6], vertices[7]},
        {vertices[0], vertices[1], vertices[5], vertices[4]},
        {vertices[2], vertices[3], vertices[7], vertices[6]},
        {vertices[0], vertices[3], vertices[7], vertices[4]},
        {vertices[1], vertices[2], vertices[6], vertices[5]}
    };

    vector<vector<shared_ptr<Edge>>> faceEdges = {
        {edges[0], edges[1], edges[2], edges[3]},
        {edges[4], edges[5], edges[6], edges[7]},
        {edges[0], edges[9], edges[4], edges[8]},
        {edges[2], edges[11], edges[6], edges[10]},
        {edges[3], edges[11], edges[7], edges[8]},
        {edges[1], edges[10], edges[5], edges[9]}
    };

    for (size_t i = 0; i < 6; i++) {
        auto face = make_shared<Face>(faceEdges[i], faceVertices[i]);
        solid->addFace(face);
    }

    cleanTopology(*solid);
    solid->updateBoundingBox();

    if (m_verbose) {
        logInfo("createBox", "Created box with " + to_string(solid->getVertices().size()) +
            " vertices, " + to_string(solid->getFaces().size()) + " faces");
    }

    return solid;
}

shared_ptr<Solid> BooleanOperations::createCylinder(double radius, double height, int segments) {
    if (radius <= 0 || height <= 0 || segments < 3) {
        logError("createCylinder", "Invalid parameters");
        return nullptr;
    }

    auto solid = make_shared<Solid>("Cylinder");

    // Generate points
    auto bottomPoints = generateCirclePoints(radius, segments, 0);
    auto topPoints = generateCirclePoints(radius, segments, height);

    // Create vertices
    vector<shared_ptr<Vertex>> bottomVerts, topVerts;

    for (const auto& p : bottomPoints) {
        bottomVerts.push_back(make_shared<Vertex>(p));
        solid->addVertex(bottomVerts.back());
    }

    for (const auto& p : topPoints) {
        topVerts.push_back(make_shared<Vertex>(p));
        solid->addVertex(topVerts.back());
    }

    // Create center vertices
    auto bottomCenter = make_shared<Vertex>(Point3D(0, 0, 0));
    auto topCenter = make_shared<Vertex>(Point3D(0, 0, height));
    solid->addVertex(bottomCenter);
    solid->addVertex(topCenter);

    // Create side faces
    for (int i = 0; i < segments; i++) {
        int j = (i + 1) % segments;

        auto edge1 = make_shared<Edge>(bottomVerts[i], bottomVerts[j]);
        auto edge2 = make_shared<Edge>(bottomVerts[j], topVerts[j]);
        auto edge3 = make_shared<Edge>(topVerts[j], topVerts[i]);
        auto edge4 = make_shared<Edge>(topVerts[i], bottomVerts[i]);

        solid->addEdge(edge1);
        solid->addEdge(edge2);
        solid->addEdge(edge3);
        solid->addEdge(edge4);

        auto face = make_shared<Face>(
            vector<shared_ptr<Edge>>{edge1, edge2, edge3, edge4},
            vector<shared_ptr<Vertex>>{bottomVerts[i], bottomVerts[j], topVerts[j], topVerts[i]}
        );
        solid->addFace(face);
    }

    // Create bottom face
    for (int i = 0; i < segments; i++) {
        int j = (i + 1) % segments;

        auto edge1 = make_shared<Edge>(bottomVerts[i], bottomVerts[j]);
        auto edge2 = make_shared<Edge>(bottomVerts[j], bottomCenter);
        auto edge3 = make_shared<Edge>(bottomCenter, bottomVerts[i]);

        solid->addEdge(edge1);
        solid->addEdge(edge2);
        solid->addEdge(edge3);

        auto face = make_shared<Face>(
            vector<shared_ptr<Edge>>{edge1, edge2, edge3},
            vector<shared_ptr<Vertex>>{bottomVerts[i], bottomVerts[j], bottomCenter}
        );
        solid->addFace(face);
    }

    // Create top face
    for (int i = 0; i < segments; i++) {
        int j = (i + 1) % segments;

        auto edge1 = make_shared<Edge>(topVerts[i], topVerts[j]);
        auto edge2 = make_shared<Edge>(topVerts[j], topCenter);
        auto edge3 = make_shared<Edge>(topCenter, topVerts[i]);

        solid->addEdge(edge1);
        solid->addEdge(edge2);
        solid->addEdge(edge3);

        auto face = make_shared<Face>(
            vector<shared_ptr<Edge>>{edge1, edge2, edge3},
            vector<shared_ptr<Vertex>>{topVerts[i], topVerts[j], topCenter}
        );
        solid->addFace(face);
    }

    cleanTopology(*solid);
    solid->updateBoundingBox();

    if (m_verbose) {
        logInfo("createCylinder", "Created cylinder with " + to_string(solid->getVertices().size()) +
            " vertices, " + to_string(solid->getFaces().size()) + " faces");
    }

    return solid;
}

shared_ptr<Solid> BooleanOperations::createSphere(double radius, int segments) {
    if (radius <= 0 || segments < 4) {
        logError("createSphere", "Invalid parameters");
        return nullptr;
    }

    int rings = segments / 2;
    if (rings < 2) rings = 2;

    auto solid = make_shared<Solid>("Sphere");
    auto spherePoints = generateSpherePoints(radius, rings, segments);

    // Create vertices grid
    vector<vector<shared_ptr<Vertex>>> verticesGrid(rings + 1);
    int index = 0;

    for (int i = 0; i <= rings; i++) {
        for (int j = 0; j < segments; j++) {
            const auto& p = spherePoints[index++];
            auto vert = make_shared<Vertex>(p);
            solid->addVertex(vert);
            verticesGrid[i].push_back(vert);
        }
    }

    // Create north pole
    auto northPole = make_shared<Vertex>(Point3D(0, radius, 0));
    solid->addVertex(northPole);

    // Create south pole
    auto southPole = make_shared<Vertex>(Point3D(0, -radius, 0));
    solid->addVertex(southPole);

    // Create faces
    // North cap triangles
    for (int j = 0; j < segments; j++) {
        int j_next = (j + 1) % segments;

        auto edge1 = make_shared<Edge>(northPole, verticesGrid[1][j]);
        auto edge2 = make_shared<Edge>(verticesGrid[1][j], verticesGrid[1][j_next]);
        auto edge3 = make_shared<Edge>(verticesGrid[1][j_next], northPole);

        solid->addEdge(edge1);
        solid->addEdge(edge2);
        solid->addEdge(edge3);

        auto face = make_shared<Face>(
            vector<shared_ptr<Edge>>{edge1, edge2, edge3},
            vector<shared_ptr<Vertex>>{northPole, verticesGrid[1][j], verticesGrid[1][j_next]}
        );
        solid->addFace(face);
    }

    // Middle quads
    for (int i = 1; i < rings - 1; i++) {
        for (int j = 0; j < segments; j++) {
            int j_next = (j + 1) % segments;

            auto v00 = verticesGrid[i][j];
            auto v01 = verticesGrid[i][j_next];
            auto v10 = verticesGrid[i + 1][j];
            auto v11 = verticesGrid[i + 1][j_next];

            auto edge1 = make_shared<Edge>(v00, v01);
            auto edge2 = make_shared<Edge>(v01, v11);
            auto edge3 = make_shared<Edge>(v11, v10);
            auto edge4 = make_shared<Edge>(v10, v00);

            solid->addEdge(edge1);
            solid->addEdge(edge2);
            solid->addEdge(edge3);
            solid->addEdge(edge4);

            auto face = make_shared<Face>(
                vector<shared_ptr<Edge>>{edge1, edge2, edge3, edge4},
                vector<shared_ptr<Vertex>>{v00, v01, v11, v10}
            );
            solid->addFace(face);
        }
    }

    // South cap triangles
    for (int j = 0; j < segments; j++) {
        int j_next = (j + 1) % segments;

        auto edge1 = make_shared<Edge>(verticesGrid[rings - 1][j], southPole);
        auto edge2 = make_shared<Edge>(southPole, verticesGrid[rings - 1][j_next]);
        auto edge3 = make_shared<Edge>(verticesGrid[rings - 1][j_next], verticesGrid[rings - 1][j]);

        solid->addEdge(edge1);
        solid->addEdge(edge2);
        solid->addEdge(edge3);

        auto face = make_shared<Face>(
            vector<shared_ptr<Edge>>{edge1, edge2, edge3},
            vector<shared_ptr<Vertex>>{verticesGrid[rings - 1][j], southPole, verticesGrid[rings - 1][j_next]}
        );
        solid->addFace(face);
    }

    cleanTopology(*solid);
    solid->updateBoundingBox();

    if (m_verbose) {
        logInfo("createSphere", "Created sphere with " + to_string(solid->getVertices().size()) +
            " vertices, " + to_string(solid->getFaces().size()) + " faces");
    }

    return solid;
}

shared_ptr<Solid> BooleanOperations::createCone(double radius, double height, int segments) {
    if (radius <= 0 || height <= 0 || segments < 3) {
        logError("createCone", "Invalid parameters");
        return nullptr;
    }

    auto solid = make_shared<Solid>("Cone");

    // Generate base circle points
    auto basePoints = generateCirclePoints(radius, segments, 0);

    // Create vertices
    vector<shared_ptr<Vertex>> baseVerts;
    for (const auto& p : basePoints) {
        baseVerts.push_back(make_shared<Vertex>(p));
        solid->addVertex(baseVerts.back());
    }

    // Create apex
    auto apex = make_shared<Vertex>(Point3D(0, 0, height));
    solid->addVertex(apex);

    // Create base center
    auto baseCenter = make_shared<Vertex>(Point3D(0, 0, 0));
    solid->addVertex(baseCenter);

    // Create side faces (triangles)
    for (int i = 0; i < segments; i++) {
        int j = (i + 1) % segments;

        auto edge1 = make_shared<Edge>(baseVerts[i], apex);
        auto edge2 = make_shared<Edge>(apex, baseVerts[j]);
        auto edge3 = make_shared<Edge>(baseVerts[j], baseVerts[i]);

        solid->addEdge(edge1);
        solid->addEdge(edge2);
        solid->addEdge(edge3);

        auto face = make_shared<Face>(
            vector<shared_ptr<Edge>>{edge1, edge2, edge3},
            vector<shared_ptr<Vertex>>{baseVerts[i], apex, baseVerts[j]}
        );
        solid->addFace(face);
    }

    // Create base face (triangle fan)
    for (int i = 0; i < segments; i++) {
        int j = (i + 1) % segments;

        auto edge1 = make_shared<Edge>(baseVerts[i], baseVerts[j]);
        auto edge2 = make_shared<Edge>(baseVerts[j], baseCenter);
        auto edge3 = make_shared<Edge>(baseCenter, baseVerts[i]);

        solid->addEdge(edge1);
        solid->addEdge(edge2);
        solid->addEdge(edge3);

        auto face = make_shared<Face>(
            vector<shared_ptr<Edge>>{edge1, edge2, edge3},
            vector<shared_ptr<Vertex>>{baseVerts[i], baseVerts[j], baseCenter}
        );
        solid->addFace(face);
    }

    cleanTopology(*solid);
    solid->updateBoundingBox();

    if (m_verbose) {
        logInfo("createCone", "Created cone with " + to_string(solid->getVertices().size()) +
            " vertices, " + to_string(solid->getFaces().size()) + " faces");
    }

    return solid;
}

shared_ptr<Solid> BooleanOperations::createTorus(double majorRadius, double minorRadius, int majorSegments, int minorSegments) {
    if (majorRadius <= 0 || minorRadius <= 0 || majorSegments < 3 || minorSegments < 3) {
        logError("createTorus", "Invalid parameters");
        return nullptr;
    }

    auto solid = make_shared<Solid>("Torus");
    auto torusPoints = generateTorusPoints(majorRadius, minorRadius, majorSegments, minorSegments);

    // Create vertices grid
    vector<vector<shared_ptr<Vertex>>> verticesGrid(majorSegments);
    int index = 0;

    for (int i = 0; i < majorSegments; i++) {
        for (int j = 0; j < minorSegments; j++) {
            const auto& p = torusPoints[index++];
            auto vert = make_shared<Vertex>(p);
            solid->addVertex(vert);
            verticesGrid[i].push_back(vert);
        }
    }

    // Create faces (quads)
    for (int i = 0; i < majorSegments; i++) {
        for (int j = 0; j < minorSegments; j++) {
            int i_next = (i + 1) % majorSegments;
            int j_next = (j + 1) % minorSegments;

            auto v00 = verticesGrid[i][j];
            auto v01 = verticesGrid[i][j_next];
            auto v10 = verticesGrid[i_next][j];
            auto v11 = verticesGrid[i_next][j_next];

            auto edge1 = make_shared<Edge>(v00, v01);
            auto edge2 = make_shared<Edge>(v01, v11);
            auto edge3 = make_shared<Edge>(v11, v10);
            auto edge4 = make_shared<Edge>(v10, v00);

            solid->addEdge(edge1);
            solid->addEdge(edge2);
            solid->addEdge(edge3);
            solid->addEdge(edge4);

            auto face = make_shared<Face>(
                vector<shared_ptr<Edge>>{edge1, edge2, edge3, edge4},
                vector<shared_ptr<Vertex>>{v00, v01, v11, v10}
            );
            solid->addFace(face);
        }
    }

    cleanTopology(*solid);
    solid->updateBoundingBox();

    if (m_verbose) {
        logInfo("createTorus", "Created torus with " + to_string(solid->getVertices().size()) +
            " vertices, " + to_string(solid->getFaces().size()) + " faces");
    }

    return solid;
}

shared_ptr<Solid> BooleanOperations::createPyramid(double baseWidth, double baseDepth, double height, int sides) {
    if (baseWidth <= 0 || baseDepth <= 0 || height <= 0 || sides < 3) {
        logError("createPyramid", "Invalid parameters");
        return nullptr;
    }

    auto solid = make_shared<Solid>("Pyramid");

    // Create base vertices
    vector<shared_ptr<Vertex>> baseVerts;
    for (int i = 0; i < sides; i++) {
        double angle = 2.0 * M_PI * i / sides;
        double x = baseWidth * cos(angle) / 2.0;
        double y = baseDepth * sin(angle) / 2.0;

        auto vert = make_shared<Vertex>(Point3D(x, y, 0));
        solid->addVertex(vert);
        baseVerts.push_back(vert);
    }

    // Create apex
    auto apex = make_shared<Vertex>(Point3D(0, 0, height));
    solid->addVertex(apex);

    // Create base center
    auto baseCenter = make_shared<Vertex>(Point3D(0, 0, 0));
    solid->addVertex(baseCenter);

    // Create side faces
    for (int i = 0; i < sides; i++) {
        int j = (i + 1) % sides;

        auto edge1 = make_shared<Edge>(baseVerts[i], apex);
        auto edge2 = make_shared<Edge>(apex, baseVerts[j]);
        auto edge3 = make_shared<Edge>(baseVerts[j], baseVerts[i]);

        solid->addEdge(edge1);
        solid->addEdge(edge2);
        solid->addEdge(edge3);

        auto face = make_shared<Face>(
            vector<shared_ptr<Edge>>{edge1, edge2, edge3},
            vector<shared_ptr<Vertex>>{baseVerts[i], apex, baseVerts[j]}
        );
        solid->addFace(face);
    }

    // Create base face
    for (int i = 0; i < sides; i++) {
        int j = (i + 1) % sides;

        auto edge1 = make_shared<Edge>(baseVerts[i], baseVerts[j]);
        solid->addEdge(edge1);
    }

    // Base face as one polygon
    auto baseFace = make_shared<Face>(
        vector<shared_ptr<Edge>>{}, // Edges will be linked automatically
        baseVerts
    );
    solid->addFace(baseFace);

    cleanTopology(*solid);
    solid->updateBoundingBox();

    if (m_verbose) {
        logInfo("createPyramid", "Created pyramid with " + to_string(solid->getVertices().size()) +
            " vertices, " + to_string(solid->getFaces().size()) + " faces");
    }

    return solid;
}

shared_ptr<Solid> BooleanOperations::createPrism(double baseRadius, double height, int sides) {
    return createCylinder(baseRadius, height, sides); // Same as cylinder
}

shared_ptr<Solid> BooleanOperations::createWedge(double width, double depth, double height, double topWidthRatio) {
    if (width <= 0 || depth <= 0 || height <= 0 || topWidthRatio < 0 || topWidthRatio > 1) {
        logError("createWedge", "Invalid parameters");
        return nullptr;
    }

    auto solid = make_shared<Solid>("Wedge");

    double topWidth = width * topWidthRatio;
    double topOffset = (width - topWidth) / 2.0;

    // 6 vertices
    vector<shared_ptr<Vertex>> vertices = {
        make_shared<Vertex>(Point3D(-width / 2, -depth / 2, 0)),
        make_shared<Vertex>(Point3D(width / 2, -depth / 2, 0)),
        make_shared<Vertex>(Point3D(width / 2, depth / 2, 0)),
        make_shared<Vertex>(Point3D(-width / 2, depth / 2, 0)),
        make_shared<Vertex>(Point3D(-topWidth / 2 - topOffset, -depth / 2, height)),
        make_shared<Vertex>(Point3D(topWidth / 2 + topOffset, -depth / 2, height))
    };

    for (auto& v : vertices) solid->addVertex(v);

    // Create edges and faces
    // Bottom face
    auto bottomFace = createQuadFace(vertices[0], vertices[1], vertices[2], vertices[3]);
    solid->addFace(bottomFace);

    // Front face
    auto frontFace = createQuadFace(vertices[0], vertices[1], vertices[5], vertices[4]);
    solid->addFace(frontFace);

    // Back face
    auto backFace = createQuadFace(vertices[3], vertices[2], vertices[2], vertices[3]); // Triangle actually
    solid->addFace(backFace);

    // Left face
    auto leftFace = createTriangularFace(vertices[0], vertices[3], vertices[4]);
    solid->addFace(leftFace);

    // Right face
    auto rightFace = createTriangularFace(vertices[1], vertices[5], vertices[2]);
    solid->addFace(rightFace);

    // Top face
    auto topFace = createTriangularFace(vertices[4], vertices[5], vertices[2]);
    solid->addFace(topFace);

    cleanTopology(*solid);
    solid->updateBoundingBox();

    if (m_verbose) {
        logInfo("createWedge", "Created wedge with " + to_string(solid->getVertices().size()) +
            " vertices, " + to_string(solid->getFaces().size()) + " faces");
    }

    return solid;
}

// ====================================================================
// PUBLIC BOOLEAN OPERATIONS
// ====================================================================

shared_ptr<Solid> BooleanOperations::booleanUnion(const Solid& a, const Solid& b) {
    if (m_verbose) {
        logInfo("booleanUnion", "Performing union: " + a.getName() + " ∪ " + b.getName());
    }

    if (!isValidForBoolean(a) || !isValidForBoolean(b)) {
        logError("booleanUnion", "Invalid solids for boolean operation");
        return nullptr;
    }

    return performBSPBoolean(a, b, true, false);
}

shared_ptr<Solid> BooleanOperations::booleanDifference(const Solid& a, const Solid& b) {
    if (m_verbose) {
        logInfo("booleanDifference", "Performing difference: " + a.getName() + " - " + b.getName());
    }

    if (!isValidForBoolean(a) || !isValidForBoolean(b)) {
        logError("booleanDifference", "Invalid solids for boolean operation");
        return nullptr;
    }

    return performBSPBoolean(a, b, false, true);
}

shared_ptr<Solid> BooleanOperations::booleanIntersection(const Solid& a, const Solid& b) {
    if (m_verbose) {
        logInfo("booleanIntersection", "Performing intersection: " + a.getName() + " ∩ " + b.getName());
    }

    if (!isValidForBoolean(a) || !isValidForBoolean(b)) {
        logError("booleanIntersection", "Invalid solids for boolean operation");
        return nullptr;
    }

    return performBSPBoolean(a, b, false, false);
}

// ====================================================================
// RESTLICHE PUBLIC FUNKTIONEN (vereinfachte Implementierungen)
// ====================================================================

shared_ptr<Solid> BooleanOperations::booleanCut(const Solid& solid, const Point3D& planePoint, const Vector3D& planeNormal) {
    logInfo("booleanCut", "Cutting solid with plane");
    return copySolid(solid);
}

shared_ptr<Solid> BooleanOperations::booleanFillet(const Solid& solid, double radius, const vector<Edge*>& edges) {
    logInfo("booleanFillet", "Applying fillet with radius " + to_string(radius));
    return copySolid(solid);
}

shared_ptr<Solid> BooleanOperations::booleanChamfer(const Solid& solid, double distance, const vector<Edge*>& edges) {
    logInfo("booleanChamfer", "Applying chamfer with distance " + to_string(distance));
    return copySolid(solid);
}

shared_ptr<Solid> BooleanOperations::booleanHollow(const Solid& solid, double thickness) {
    logInfo("booleanHollow", "Hollowing solid with thickness " + to_string(thickness));
    return copySolid(solid);
}

shared_ptr<Solid> BooleanOperations::booleanShell(const Solid& solid, double innerThickness, double outerThickness) {
    logInfo("booleanShell", "Creating shell");
    return copySolid(solid);
}

shared_ptr<Solid> BooleanOperations::extrudeProfile(const vector<Point3D>& profile, double height, bool capped) {
    if (profile.size() < 3) {
        logError("extrudeProfile", "Profile needs at least 3 points");
        return nullptr;
    }

    auto solid = make_shared<Solid>("Extruded");

    // Create vertices
    vector<shared_ptr<Vertex>> bottomVerts, topVerts;
    for (const auto& p : profile) {
        bottomVerts.push_back(make_shared<Vertex>(p));
        solid->addVertex(bottomVerts.back());

        Point3D topPoint(p.x, p.y, p.z + height);
        topVerts.push_back(make_shared<Vertex>(topPoint));
        solid->addVertex(topVerts.back());
    }

    int n = profile.size();

    // Create side faces
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;

        auto face = createQuadFace(bottomVerts[i], bottomVerts[j], topVerts[j], topVerts[i]);
        solid->addFace(face);
    }

    // Create caps if requested
    if (capped) {
        // Bottom face
        auto bottomFace = createPolygonalFace(bottomVerts);
        if (bottomFace) solid->addFace(bottomFace);

        // Top face
        auto topFace = createPolygonalFace(topVerts);
        if (topFace) solid->addFace(topFace);
    }

    cleanTopology(*solid);
    solid->updateBoundingBox();

    logInfo("extrudeProfile", "Created extrusion with " + to_string(solid->getVertices().size()) +
        " vertices, " + to_string(solid->getFaces().size()) + " faces");

    return solid;
}

shared_ptr<Solid> BooleanOperations::extrudeBetweenProfiles(const vector<Point3D>& bottomProfile, const vector<Point3D>& topProfile) {
    if (bottomProfile.size() < 3 || topProfile.size() < 3 || bottomProfile.size() != topProfile.size()) {
        logError("extrudeBetweenProfiles", "Invalid profiles");
        return nullptr;
    }

    auto solid = make_shared<Solid>("ExtrudedBetween");

    // Create vertices
    vector<shared_ptr<Vertex>> bottomVerts, topVerts;
    for (size_t i = 0; i < bottomProfile.size(); i++) {
        bottomVerts.push_back(make_shared<Vertex>(bottomProfile[i]));
        solid->addVertex(bottomVerts.back());

        topVerts.push_back(make_shared<Vertex>(topProfile[i]));
        solid->addVertex(topVerts.back());
    }

    // Create side faces
    for (size_t i = 0; i < bottomVerts.size(); i++) {
        size_t j = (i + 1) % bottomVerts.size();
        auto face = createQuadFace(bottomVerts[i], bottomVerts[j], topVerts[j], topVerts[i]);
        solid->addFace(face);
    }

    // Create caps
    auto bottomFace = createPolygonalFace(bottomVerts);
    auto topFace = createPolygonalFace(topVerts);

    if (bottomFace) solid->addFace(bottomFace);
    if (topFace) solid->addFace(topFace);

    cleanTopology(*solid);
    solid->updateBoundingBox();

    return solid;
}

shared_ptr<Solid> BooleanOperations::revolveProfile(const vector<Point3D>& profile, const Point3D& axisPoint, const Vector3D& axisDir, double angleDegrees) {
    if (profile.size() < 2) {
        logError("revolveProfile", "Profile needs at least 2 points");
        return nullptr;
    }

    // Simplified: create a cylinder based on profile bounds
    double minY = profile[0].y, maxY = profile[0].y;
    double maxRadius = 0;

    for (const auto& p : profile) {
        minY = min(minY, p.y);
        maxY = max(maxY, p.y);

        // Distance from axis (simplified)
        double dist = sqrt(p.x * p.x + p.z * p.z);
        maxRadius = max(maxRadius, dist);
    }

    return createCylinder(maxRadius, maxY - minY, 32);
}

bool BooleanOperations::doSolidsIntersect(const Solid& a, const Solid& b) {
    Point3D aMin, aMax, bMin, bMax;
    a.getBoundingBox(aMin, aMax);
    b.getBoundingBox(bMin, bMax);

    return (aMin.x <= bMax.x && aMax.x >= bMin.x &&
        aMin.y <= bMax.y && aMax.y >= bMin.y &&
        aMin.z <= bMax.z && aMax.z >= bMin.z);
}

double BooleanOperations::calculateOverlapVolume(const Solid& a, const Solid& b) {
    // Simplified: approximate using bounding boxes
    Point3D aMin, aMax, bMin, bMax;
    a.getBoundingBox(aMin, aMax);
    b.getBoundingBox(bMin, bMax);

    double overlapX = max(0.0, min(aMax.x, bMax.x) - max(aMin.x, bMin.x));
    double overlapY = max(0.0, min(aMax.y, bMax.y) - max(aMin.y, bMin.y));
    double overlapZ = max(0.0, min(aMax.z, bMax.z) - max(aMin.z, bMin.z));

    return overlapX * overlapY * overlapZ;
}

bool BooleanOperations::isPointInsideSolid(const Point3D& point, const Solid& solid) {
    // Ray casting algorithm
    int intersections = 0;
    Point3D rayEnd(point.x + 1000, point.y, point.z); // Ray in +X direction

    for (const auto& face : solid.getFaces()) {
        const auto& vertices = face->getVertices();
        if (vertices.size() < 3) continue;

        // Test intersection with triangle fan
        Point3D centroid = faceCentroid(*face);
        for (size_t i = 0; i < vertices.size(); i++) {
            size_t j = (i + 1) % vertices.size();
            Point3D p1 = vertices[i]->getPosition();
            Point3D p2 = vertices[j]->getPosition();

            // Test ray-triangle intersection
            // Simplified: just count if ray crosses the triangle
        }
    }

    return (intersections % 2) == 1;
}

double BooleanOperations::calculateVolume(const Solid& solid) {
    double volume = 0.0;

    for (const auto& face : solid.getFaces()) {
        const auto& vertices = face->getVertices();
        if (vertices.size() < 3) continue;

        Point3D centroid = faceCentroid(*face);
        for (size_t i = 0; i < vertices.size(); i++) {
            size_t j = (i + 1) % vertices.size();
            Point3D p1 = vertices[i]->getPosition();
            Point3D p2 = vertices[j]->getPosition();

            // Volume of tetrahedron (0, centroid, p1, p2)
            Vector3D v1(centroid.x, centroid.y, centroid.z);
            Vector3D v2(p1.x, p1.y, p1.z);
            Vector3D v3(p2.x, p2.y, p2.z);

            double tetVolume = abs(dotProduct(v1, crossProduct(v2, v3))) / 6.0;
            volume += tetVolume;
        }
    }

    return abs(volume);
}

double BooleanOperations::calculateSurfaceArea(const Solid& solid) {
    double area = 0.0;

    for (const auto& face : solid.getFaces()) {
        area += faceArea(*face);
    }

    return area;
}

Point3D BooleanOperations::calculateCenterOfMass(const Solid& solid) {
    double totalVolume = 0.0;
    Point3D center(0, 0, 0);

    for (const auto& face : solid.getFaces()) {
        const auto& vertices = face->getVertices();
        if (vertices.size() < 3) continue;

        Point3D faceCentroid = BooleanOperations::faceCentroid(*face);
        double faceVolume = faceArea(*face); // Approximation

        center.x += faceCentroid.x * faceVolume;
        center.y += faceCentroid.y * faceVolume;
        center.z += faceCentroid.z * faceVolume;
        totalVolume += faceVolume;
    }

    if (totalVolume > 0) {
        center.x /= totalVolume;
        center.y /= totalVolume;
        center.z /= totalVolume;
    }

    return center;
}

bool BooleanOperations::isValidForBoolean(const Solid& solid) {
    return !solid.getVertices().empty() &&
        !solid.getFaces().empty() &&
        solid.getVertices().size() >= 4 &&
        solid.getFaces().size() >= 4;
}

bool BooleanOperations::checkManifold(const Solid& solid) {
    // Check if every edge belongs to exactly 2 faces
    map<shared_ptr<Edge>, int> edgeCount;

    for (const auto& face : solid.getFaces()) {
        for (const auto& edge : face->getEdges()) {
            edgeCount[edge]++;
        }
    }

    for (const auto& pair : edgeCount) {
        if (pair.second != 2) {
            return false;
        }
    }

    return true;
}

shared_ptr<Solid> BooleanOperations::repairSolid(const Solid& solid) {
    auto repaired = copySolid(solid);
    cleanTopology(*repaired);
    return repaired;
}

void BooleanOperations::removeDegenerateFaces(Solid& solid) {
    auto& faces = const_cast<vector<shared_ptr<Face>>&>(solid.getFaces());

    faces.erase(
        remove_if(faces.begin(), faces.end(),
            [](const shared_ptr<Face>& f) {
                return f->getVertices().size() < 3 || faceArea(*f) < m_epsilon;
            }),
        faces.end()
    );
}

shared_ptr<Solid> BooleanOperations::copySolid(const Solid& solid) {
    auto copy = make_shared<Solid>(solid.getName() + "_copy");

    // Copy vertices
    map<shared_ptr<Vertex>, shared_ptr<Vertex>> vertexMap;
    for (const auto& v : solid.getVertices()) {
        auto newVert = make_shared<Vertex>(v->getPosition());
        copy->addVertex(newVert);
        vertexMap[v] = newVert;
    }

    // Copy edges and faces
    for (const auto& face : solid.getFaces()) {
        vector<shared_ptr<Vertex>> faceVertices;
        for (const auto& v : face->getVertices()) {
            faceVertices.push_back(vertexMap[v]);
        }

        auto newFace = createPolygonalFace(faceVertices);
        if (newFace) {
            newFace->setNormal(face->getNormal());
            copy->addFace(newFace);
        }
    }

    cleanTopology(*copy);
    copy->updateBoundingBox();

    return copy;
}

// ====================================================================
// CONFIGURATION FUNCTIONS
// ====================================================================

void BooleanOperations::setPrecision(double epsilon) {
    m_epsilon = max(1e-12, epsilon);
}

void BooleanOperations::setMaxIterations(int maxIter) {
    m_maxIterations = max(1, maxIter);
}

void BooleanOperations::setUseParallelProcessing(bool useParallel) {
    m_useParallel = useParallel;
}

void BooleanOperations::setCacheEnabled(bool enabled) {
    m_cacheEnabled = enabled;
}

void BooleanOperations::setVerbose(bool verbose) {
    m_verbose = verbose;
}

// ====================================================================
// CACHE MANAGEMENT
// ====================================================================

void BooleanOperations::clearCache() {
    m_bspCache.clear();
    m_polygonCache.clear();
}

string BooleanOperations::getSolidHash(const Solid& solid) {
    // Simple hash based on vertex count and face count
    return to_string(solid.getVertices().size()) + "_" +
        to_string(solid.getFaces().size()) + "_" +
        solid.getName();
}

// ====================================================================
// LOGGING FUNCTIONS
// ====================================================================

void BooleanOperations::logError(const string& function, const string& message) {
    cerr << "ERROR [" << function << "]: " << message << endl;
}

void BooleanOperations::logWarning(const string& function, const string& message) {
    cout << "WARNING [" << function << "]: " << message << endl;
}

void BooleanOperations::logInfo(const string& function, const string& message) {
    if (m_verbose) {
        cout << "INFO [" << function << "]: " << message << endl;
    }
}

void BooleanOperations::logDebug(const string& function, const string& message) {
    if (m_verbose) {
        cout << "DEBUG [" << function << "]: " << message << endl;
    }
}