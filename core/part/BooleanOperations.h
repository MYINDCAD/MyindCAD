#pragma once

#include <memory>
#include <vector>
#include <string>

// Forward declarations
class Solid;
class Face;
class Edge;
class Vertex;
struct Point3D;
struct Vector3D;

class BooleanOperations {
public:
    // ===== BOOLEAN OPERATIONS (CSG) =====
    static std::shared_ptr<Solid> booleanUnion(const Solid& a, const Solid& b);
    static std::shared_ptr<Solid> booleanDifference(const Solid& a, const Solid& b);
    static std::shared_ptr<Solid> booleanIntersection(const Solid& a, const Solid& b);

    // Advanced boolean operations
    static std::shared_ptr<Solid> booleanCut(const Solid& solid, const Point3D& planePoint,
        const Vector3D& planeNormal);
    static std::shared_ptr<Solid> booleanFillet(const Solid& solid, double radius);
    static std::shared_ptr<Solid> booleanChamfer(const Solid& solid, double distance);

    // ===== EXTRUSION OPERATIONS =====
    static std::shared_ptr<Solid> extrudeProfile(const std::vector<Point3D>& profile,
        double height, bool capped = true);
    static std::shared_ptr<Solid> extrudeBetweenProfiles(const std::vector<Point3D>& bottomProfile,
        const std::vector<Point3D>& topProfile);

    // ===== REVOLUTION OPERATIONS =====
    static std::shared_ptr<Solid> revolveProfile(const std::vector<Point3D>& profile,
        const Point3D& axisPoint,
        const Vector3D& axisDir,
        double angleDegrees = 360.0);
    static std::shared_ptr<Solid> revolveAroundEdge(const std::vector<Point3D>& profile,
        const Edge& axisEdge,
        double angleDegrees = 360.0);

    // ===== GEOMETRY ANALYSIS =====
    static bool doSolidsIntersect(const Solid& a, const Solid& b);
    static double calculateOverlapVolume(const Solid& a, const Solid& b);
    static std::vector<Point3D> findFaceIntersection(const Face& faceA, const Face& faceB);
    static bool isPointInsideSolid(const Point3D& point, const Solid& solid);

    // ===== VALIDATION & REPAIR =====
    static bool isValidForBoolean(const Solid& solid);
    static bool checkManifold(const Solid& solid);
    static std::shared_ptr<Solid> repairSolid(const Solid& solid);
    static void removeDegenerateFaces(Solid& solid);

    // ===== CONSTRAINT SOLVER INTEGRATION =====
    static std::shared_ptr<Solid> applyConstraints(Solid& solid,
        const std::vector<std::string>& constraints);
    static bool checkGeometricConstraints(const Solid& solid,
        const std::vector<std::string>& constraints);

    // ===== UTILITY FUNCTIONS =====
    static std::shared_ptr<Solid> mergeSolids(const std::vector<std::shared_ptr<Solid>>& solids);
    static std::shared_ptr<Solid> copySolid(const Solid& solid);
    static void transformSolid(Solid& solid, const std::vector<double>& transformationMatrix);

    // ===== DIAGNOSTICS =====
    static void printSolidInfo(const Solid& solid);
    static bool validateTopology(const Solid& solid);
    static std::string getSolidStats(const Solid& solid);

    // Prevent instantiation
    BooleanOperations() = delete;
    ~BooleanOperations() = delete;

private:
    // ===== PRIVATE HELPER FUNCTIONS =====

    // Boolean operation helpers
    static std::vector<Point3D> clipFaceAgainstPlane(const Face& face,
        const Point3D& planePoint,
        const Vector3D& planeNormal);
    static std::shared_ptr<Face> createFaceFromPoints(const std::vector<Point3D>& points);
    static std::vector<std::shared_ptr<Edge>> findIntersectionEdges(const Solid& a, const Solid& b);

    // Extrusion helpers
    static std::vector<Point3D> offsetProfile(const std::vector<Point3D>& profile,
        double offset);
    static std::shared_ptr<Face> createExtrusionFace(const std::vector<Point3D>& bottom,
        const std::vector<Point3D>& top,
        int startIdx, int endIdx);

    // Revolution helpers
    static Point3D rotatePoint(const Point3D& point, const Point3D& axisPoint,
        const Vector3D& axisDir, double angleDegrees);
    static std::shared_ptr<Face> createRevolvedFace(const std::vector<Point3D>& profile,
        int profileIndex,
        const Point3D& axisPoint,
        const Vector3D& axisDir,
        double angleStep);

    // Geometry helpers
    static double pointToPlaneDistance(const Point3D& point,
        const Point3D& planePoint,
        const Vector3D& planeNormal);
    static Vector3D faceNormal(const Face& face);
    static double faceArea(const Face& face);
    static Point3D faceCentroid(const Face& face);

    // Constraint solver helpers
    static bool applyDistanceConstraint(Solid& solid, const std::string& constraint);
    static bool applyAngleConstraint(Solid& solid, const std::string& constraint);
    static bool applyParallelConstraint(Solid& solid, const std::string& constraint);
    static bool applyPerpendicularConstraint(Solid& solid, const std::string& constraint);

    // Topology helpers
    static void rebuildFaceEdges(Face& face);
    static void fixVertexConnectivity(Solid& solid);
    static void removeDuplicateVertices(Solid& solid);
    static void orientFacesConsistently(Solid& solid);

    // Debug helpers
    static void logIntersection(const std::string& operation,
        const Solid& a, const Solid& b,
        const std::vector<Point3D>& intersectionPoints);
};