#pragma once

#include <memory>
#include <vector>
#include <string>
#include <map>
#include <set>

class Solid;
class Face;
class Edge;
class Vertex;
struct Point3D;
struct Vector3D;

class BooleanOperations {
public:
    // ===== PRIMITIVE CREATION =====
    static std::shared_ptr<Solid> createBox(double width, double height, double depth);
    static std::shared_ptr<Solid> createCylinder(double radius, double height, int segments = 32);
    static std::shared_ptr<Solid> createSphere(double radius, int segments = 24);
    static std::shared_ptr<Solid> createCone(double radius, double height, int segments = 32);
    static std::shared_ptr<Solid> createTorus(double majorRadius, double minorRadius, int majorSegments = 32, int minorSegments = 16);
    static std::shared_ptr<Solid> createPyramid(double baseWidth, double baseDepth, double height, int sides = 4);
    static std::shared_ptr<Solid> createPrism(double baseRadius, double height, int sides = 6);
    static std::shared_ptr<Solid> createWedge(double width, double depth, double height, double topWidthRatio = 0.5);

    // ===== ADVANCED PRIMITIVES =====
    static std::shared_ptr<Solid> createHelix(double radius, double pitch, double height, int turns = 3, int segmentsPerTurn = 12);
    static std::shared_ptr<Solid> createSpiral(double startRadius, double endRadius, double height, int turns = 5, int segmentsPerTurn = 16);
    static std::shared_ptr<Solid> createGear(double pitchRadius, double toothHeight, int teeth, double thickness);
    static std::shared_ptr<Solid> createSpring(double coilRadius, double wireRadius, double height, int coils = 10, int segmentsPerCoil = 12);

    // ===== BOOLEAN OPERATIONS =====
    static std::shared_ptr<Solid> booleanUnion(const Solid& a, const Solid& b);
    static std::shared_ptr<Solid> booleanDifference(const Solid& a, const Solid& b);
    static std::shared_ptr<Solid> booleanIntersection(const Solid& a, const Solid& b);

    static std::shared_ptr<Solid> booleanCut(const Solid& solid, const Point3D& planePoint, const Vector3D& planeNormal);
    static std::shared_ptr<Solid> booleanFillet(const Solid& solid, double radius, const std::vector<Edge*>& edges = {});
    static std::shared_ptr<Solid> booleanChamfer(const Solid& solid, double distance, const std::vector<Edge*>& edges = {});
    static std::shared_ptr<Solid> booleanHollow(const Solid& solid, double thickness);
    static std::shared_ptr<Solid> booleanShell(const Solid& solid, double innerThickness, double outerThickness = 0);

    // ===== EXTRUSION OPERATIONS =====
    static std::shared_ptr<Solid> extrudeProfile(const std::vector<Point3D>& profile, double height, bool capped = true);
    static std::shared_ptr<Solid> extrudeBetweenProfiles(const std::vector<Point3D>& bottomProfile, const std::vector<Point3D>& topProfile);
    static std::shared_ptr<Solid> extrudeAlongPath(const std::vector<Point3D>& profile, const std::vector<Point3D>& path, bool capped = true);
    static std::shared_ptr<Solid> extrudeWithTaper(const std::vector<Point3D>& profile, double height, double taperAngle);
    static std::shared_ptr<Solid> extrudeWithRotation(const std::vector<Point3D>& profile, double height, double twistAngle);

    // ===== REVOLUTION OPERATIONS =====
    static std::shared_ptr<Solid> revolveProfile(const std::vector<Point3D>& profile, const Point3D& axisPoint, const Vector3D& axisDir, double angleDegrees = 360.0);
    static std::shared_ptr<Solid> revolveAroundEdge(const std::vector<Point3D>& profile, const Edge& axisEdge, double angleDegrees = 360.0);
    static std::shared_ptr<Solid> revolveWithTwist(const std::vector<Point3D>& profile, const Point3D& axisPoint, const Vector3D& axisDir, double angleDegrees, double twistAngle);

    // ===== SWEEP & LOFT OPERATIONS =====
    static std::shared_ptr<Solid> createLoft(const std::vector<std::vector<Point3D>>& profiles, bool ruled = true, bool smooth = false);
    static std::shared_ptr<Solid> createSweep(const std::vector<Point3D>& profile, const std::vector<Point3D>& path, bool frenet = true, bool closed = false);
    static std::shared_ptr<Solid> createCoonsPatch(const std::vector<Point3D>& boundaryCurves);
    static std::shared_ptr<Solid> createRuledSurface(const std::vector<Point3D>& curve1, const std::vector<Point3D>& curve2);

    // ===== TRANSFORM OPERATIONS =====
    static std::shared_ptr<Solid> scaleSolid(const Solid& solid, double scaleX, double scaleY, double scaleZ);
    static std::shared_ptr<Solid> translateSolid(const Solid& solid, double dx, double dy, double dz);
    static std::shared_ptr<Solid> rotateSolid(const Solid& solid, const Point3D& center, const Vector3D& axis, double angleDegrees);
    static std::shared_ptr<Solid> mirrorSolid(const Solid& solid, const Point3D& planePoint, const Vector3D& planeNormal);
    static std::shared_ptr<Solid> arraySolid(const Solid& solid, int countX, int countY, int countZ, double spacingX, double spacingY, double spacingZ);
    static std::shared_ptr<Solid> patternSolid(const Solid& solid, const std::vector<Point3D>& positions);

    // ===== COMBINATION OPERATIONS =====
    static std::shared_ptr<Solid> mergeSolids(const std::vector<std::shared_ptr<Solid>>& solids);
    static std::shared_ptr<Solid> splitSolid(const Solid& solid, const Point3D& planePoint, const Vector3D& planeNormal);
    static std::shared_ptr<Solid> assembleSolids(const std::vector<std::shared_ptr<Solid>>& solids, const std::vector<std::string>& constraints);

    // ===== GEOMETRY ANALYSIS =====
    static bool doSolidsIntersect(const Solid& a, const Solid& b);
    static double calculateOverlapVolume(const Solid& a, const Solid& b);
    static std::vector<Point3D> findFaceIntersection(const Face& faceA, const Face& faceB);
    static bool isPointInsideSolid(const Point3D& point, const Solid& solid);
    static double calculateVolume(const Solid& solid);
    static double calculateSurfaceArea(const Solid& solid);
    static Point3D calculateCenterOfMass(const Solid& solid);
    static Vector3D calculatePrincipalAxes(const Solid& solid);
    static double calculateMass(const Solid& solid, double density = 1.0);
    static std::vector<double> calculateMomentsOfInertia(const Solid& solid, double density = 1.0);

    // ===== TOPOLOGY ANALYSIS =====
    static bool isValidForBoolean(const Solid& solid);
    static bool checkManifold(const Solid& solid);
    static bool isWatertight(const Solid& solid);
    static int countHoles(const Solid& solid);
    static int countShells(const Solid& solid);
    static std::vector<std::shared_ptr<Solid>> splitIntoConnectedComponents(const Solid& solid);

    // ===== VALIDATION & REPAIR =====
    static std::shared_ptr<Solid> repairSolid(const Solid& solid);
    static void removeDegenerateFaces(Solid& solid);
    static void fixNonManifoldEdges(Solid& solid);
    static void removeSelfIntersections(Solid& solid);
    static void closeHoles(Solid& solid);
    static void fixNormals(Solid& solid);
    static void unifyNormals(Solid& solid);
    static std::shared_ptr<Solid> makeManifold(const Solid& solid);

    // ===== MESH OPERATIONS =====
    static std::shared_ptr<Solid> remeshSolid(const Solid& solid, double targetEdgeLength);
    static std::shared_ptr<Solid> smoothSolid(const Solid& solid, int iterations = 5, double lambda = 0.5);
    static std::shared_ptr<Solid> decimateSolid(const Solid& solid, double reductionFactor);
    static std::shared_ptr<Solid> refineSolid(const Solid& solid, int subdivisionLevel = 1);
    static std::shared_ptr<Solid> subdivideSolid(const Solid& solid, const std::string& scheme = "loop");
    static std::shared_ptr<Solid> simplifySolid(const Solid& solid, double tolerance);

    // ===== FEATURE RECOGNITION =====
    static std::vector<Edge*> findSharpEdges(const Solid& solid, double angleThreshold = 30.0);
    static std::vector<Face*> findPlanarFaces(const Solid& solid, double tolerance = 1e-3);
    static std::vector<Face*> findCylindricalFaces(const Solid& solid, double tolerance = 1e-3);
    static std::vector<Face*> findSphericalFaces(const Solid& solid, double tolerance = 1e-3);
    static std::vector<Face*> findConicalFaces(const Solid& solid, double tolerance = 1e-3);
    static std::vector<Face*> findToroidalFaces(const Solid& solid, double tolerance = 1e-3);

    // ===== CONSTRAINT OPERATIONS =====
    static std::shared_ptr<Solid> applyConstraints(Solid& solid, const std::vector<std::string>& constraints);
    static bool checkGeometricConstraints(const Solid& solid, const std::vector<std::string>& constraints);
    static std::vector<std::string> findViolatedConstraints(const Solid& solid, const std::vector<std::string>& constraints);
    static std::shared_ptr<Solid> solveConstraints(Solid& solid, const std::vector<std::string>& constraints);

    // ===== UTILITY FUNCTIONS =====
    static std::shared_ptr<Solid> copySolid(const Solid& solid);
    static std::shared_ptr<Solid> deepCopySolid(const Solid& solid);
    static void transformSolid(Solid& solid, const std::vector<double>& transformationMatrix);
    static std::shared_ptr<Solid> makeThickness(const Solid& solid, double thickness);
    static std::shared_ptr<Solid> offsetSolid(const Solid& solid, double offset);
    static std::shared_ptr<Solid> inflateSolid(const Solid& solid, double distance);
    static std::shared_ptr<Solid> deflateSolid(const Solid& solid, double distance);

    // ===== IMPORT/EXPORT HELPERS =====
    static std::shared_ptr<Solid> importFromFile(const std::string& filename, const std::string& format = "auto");
    static bool exportToFile(const Solid& solid, const std::string& filename, const std::string& format = "STL");
    static std::string getSolidInfo(const Solid& solid, bool detailed = false);
    static std::vector<std::string> getSupportedFormats();

    // ===== DIAGNOSTICS =====
    static void printSolidInfo(const Solid& solid);
    static bool validateTopology(const Solid& solid);
    static std::string getSolidStats(const Solid& solid);
    static void visualizeSolid(const Solid& solid);
    static void debugSolid(const Solid& solid, const std::string& filename = "debug_output.txt");
    static std::map<std::string, std::string> getDiagnostics(const Solid& solid);

    // ===== PERFORMANCE OPTIONS =====
    static void setPrecision(double epsilon);
    static void setMaxIterations(int maxIter);
    static void setUseParallelProcessing(bool useParallel);
    static void setCacheEnabled(bool enabled);
    static void setVerbose(bool verbose);

    // Prevent instantiation
    BooleanOperations() = delete;
    ~BooleanOperations() = delete;

    // ===== CONSTANTS =====
    static constexpr double DEFAULT_EPSILON = 1e-6;
    static constexpr int DEFAULT_SEGMENTS = 32;
    static constexpr int MAX_ITERATIONS = 1000;

private:
    // ===== PRIVATE STRUCTURES =====
    struct Plane {
        Point3D point;
        Vector3D normal;
        double constant;

        Plane(const Point3D& p, const Vector3D& n);
        double signedDistance(const Point3D& p) const;
        bool isPointOnPlane(const Point3D& p, double epsilon = DEFAULT_EPSILON) const;
        Point3D projectPoint(const Point3D& p) const;
    };

    struct Polygon {
        std::vector<Point3D> vertices;
        Plane plane;
        std::vector<int> edges;

        Polygon(const std::vector<Point3D>& verts);
        Polygon(const std::vector<Point3D>& verts, const Plane& pl);
        bool isValid() const;
        bool isConvex() const;
        double area() const;
        Point3D centroid() const;
        std::vector<Polygon> triangulate() const;
        std::vector<Point3D> getPoints() const { return vertices; }
    };

    struct BSPNode {
        Plane plane;
        std::vector<Polygon> polygons;
        std::unique_ptr<BSPNode> front;
        std::unique_ptr<BSPNode> back;

        BSPNode(const std::vector<Polygon>& polyList);
        std::vector<Polygon> clipPolygons(const std::vector<Polygon>& polygons) const;
        void invert();
        bool isConvex() const;
        void getAllPolygons(std::vector<Polygon>& result) const;
    };

    struct IntersectionInfo {
        Point3D point;
        double parameter;
        bool valid;
        Edge* edge1;
        Edge* edge2;
        Face* face1;
        Face* face2;
    };

    // ===== PRIVATE HELPER FUNCTIONS =====

    // Boolean operation helpers
    static std::vector<Point3D> clipPolygon(const std::vector<Point3D>& polygon, const Plane& clipPlane);
    static std::vector<Point3D> clipFaceAgainstPlane(const Face& face, const Plane& plane);
    static std::shared_ptr<Face> createFaceFromPoints(const std::vector<Point3D>& points);
    static std::vector<std::shared_ptr<Edge>> findIntersectionEdges(const Solid& a, const Solid& b);
    static std::shared_ptr<BSPNode> buildBSPTree(const Solid& solid);
    static std::vector<Polygon> solidToPolygons(const Solid& solid);
    static std::shared_ptr<Solid> polygonsToSolid(const std::vector<Polygon>& polygons, const std::string& name);
    static std::shared_ptr<Solid> performBSPBoolean(const Solid& a, const Solid& b, bool isUnion, bool isDifference);

    // Primitive helpers
    static std::vector<Point3D> generateCirclePoints(double radius, int segments, double z = 0.0, double startAngle = 0.0);
    static std::vector<Point3D> generateArcPoints(double radius, double startAngle, double endAngle, int segments, double z = 0.0);
    static std::vector<Point3D> generateSpherePoints(double radius, int rings, int segments);
    static std::vector<Point3D> generateTorusPoints(double majorRadius, double minorRadius, int majorSegments, int minorSegments);
    static std::shared_ptr<Face> createPolygonalFace(const std::vector<std::shared_ptr<Vertex>>& vertices);
    static std::shared_ptr<Face> createTriangularFace(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3);
    static std::shared_ptr<Face> createQuadFace(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3, std::shared_ptr<Vertex> v4);

    // Geometry helpers
    static double pointToPlaneDistance(const Point3D& point, const Point3D& planePoint, const Vector3D& planeNormal);
    static Vector3D faceNormal(const Face& face);
    static double faceArea(const Face& face);
    static Point3D faceCentroid(const Face& face);
    static bool pointInPolygon(const Point3D& point, const std::vector<Point3D>& polygon, const Vector3D& normal);
    static bool linesIntersect(const Point3D& p1, const Point3D& p2, const Point3D& q1, const Point3D& q2, Point3D& intersection);
    static double triangleArea(const Point3D& a, const Point3D& b, const Point3D& c);
    static double polygonArea(const std::vector<Point3D>& polygon);
    static bool isConvexPolygon(const std::vector<Point3D>& polygon);
    static std::vector<Point3D> convexHull(const std::vector<Point3D>& points);

    // Constraint helpers
    static bool applyDistanceConstraint(Solid& solid, const std::string& constraint);
    static bool applyAngleConstraint(Solid& solid, const std::string& constraint);
    static bool applyParallelConstraint(Solid& solid, const std::string& constraint);
    static bool applyPerpendicularConstraint(Solid& solid, const std::string& constraint);
    static bool applyConcentricConstraint(Solid& solid, const std::string& constraint);
    static bool applyTangentConstraint(Solid& solid, const std::string& constraint);
    static bool applyCoincidentConstraint(Solid& solid, const std::string& constraint);
    static bool applyMidpointConstraint(Solid& solid, const std::string& constraint);

    // Topology helpers
    static void rebuildFaceEdges(Face& face);
    static void fixVertexConnectivity(Solid& solid);
    static void removeDuplicateVertices(Solid& solid);
    static void orientFacesConsistently(Solid& solid);
    static void computeFaceNormals(Solid& solid);
    static void ensureWindingOrder(Face& face);
    static std::vector<std::shared_ptr<Vertex>> getFaceVerticesInOrder(const Face& face);
    static void linkEdgesToFaces(Solid& solid);
    static void linkVerticesToEdges(Solid& solid);
    static void cleanTopology(Solid& solid);

    // Mesh operation helpers
    static std::vector<Point3D> laplacianSmooth(const std::vector<Point3D>& vertices, const std::vector<std::pair<int, int>>& edges, double lambda = 0.5);
    static std::vector<std::pair<int, int>> getEdgeList(const Solid& solid);
    static std::vector<std::vector<int>> getFaceVertexIndices(const Solid& solid);
    static std::shared_ptr<Solid> subdivideCatmullClark(const Solid& solid);
    static std::shared_ptr<Solid> subdivideLoop(const Solid& solid);

    // Feature recognition helpers
    static double edgeDihedralAngle(const Edge& edge);
    static bool isEdgeSharp(const Edge& edge, double threshold);
    static bool isFacePlanar(const Face& face, double tolerance);
    static bool isFaceCylindrical(const Face& face, double tolerance, double& radius, Vector3D& axis);
    static bool isFaceSpherical(const Face& face, double tolerance, double& radius, Point3D& center);

    // Math helpers
    static Point3D rotatePoint(const Point3D& point, const Point3D& center, const Vector3D& axis, double angleDegrees);
    static Vector3D crossProduct(const Vector3D& a, const Vector3D& b);
    static double dotProduct(const Vector3D& a, const Vector3D& b);
    static double vectorLength(const Vector3D& v);
    static Vector3D normalize(const Vector3D& v);
    static double angleBetweenVectors(const Vector3D& a, const Vector3D& b);
    static Point3D linePlaneIntersection(const Point3D& linePoint, const Vector3D& lineDir, const Plane& plane);

    // Debug helpers
    static void logError(const std::string& function, const std::string& message);
    static void logWarning(const std::string& function, const std::string& message);
    static void logInfo(const std::string& function, const std::string& message);
    static void logDebug(const std::string& function, const std::string& message);

    // ===== STATIC CONFIGURATION =====
    static double m_epsilon;
    static int m_maxIterations;
    static bool m_useParallel;
    static bool m_cacheEnabled;
    static bool m_verbose;

    // ===== CACHE =====
    static std::map<std::string, std::shared_ptr<BSPNode>> m_bspCache;
    static std::map<std::string, std::vector<Polygon>> m_polygonCache;
    static void clearCache();
    static std::string getSolidHash(const Solid& solid);
};