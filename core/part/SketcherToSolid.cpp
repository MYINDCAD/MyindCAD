// SketcherToSolid.cpp
#include "SketcherToSolid.h"
#include "solid.h"
#include "Face.h"
#include "Edge.h"
#include "Vertex.h"
#include "Point3D.h"
#include "Sketch.h"
#include "curve.h"
#include <cmath>
#include <vector>
#include <stdexcept>

using namespace std;

Solid* SketcherToSolid::extrudeSketch(const Sketch* sketch, double height, bool bothSides) {
    if (!sketch || height <= 0) {
        throw invalid_argument("Invalid sketch or height");
    }

    // Geschlossenes Profil vom Sketch erhalten
    vector<Point2D> profile = sketch->getClosedProfile();
    if (profile.size() < 3) {
        throw runtime_error("Sketch does not have a closed profile");
    }

    Solid* solid = new Solid();

    // Vertices erstellen (unten und oben)
    vector<Vertex*> bottomVertices;
    vector<Vertex*> topVertices;

    double zBottom = bothSides ? -height / 2 : 0;
    double zTop = bothSides ? height / 2 : height;

    for (const Point2D& p : profile) {
        Vertex* bottom = new Vertex(Point3D(p.x, p.y, zBottom));
        Vertex* top = new Vertex(Point3D(p.x, p.y, zTop));

        bottomVertices.push_back(bottom);
        topVertices.push_back(top);

        solid->addVertex(bottom);
        solid->addVertex(top);
    }

    int n = profile.size();

    // Seitenflächen (Quads) erstellen
    for (int i = 0; i < n; i++) {
        int next = (i + 1) % n;

        // Vertikale Kanten
        Edge* sideEdge1 = new Edge(bottomVertices[i], topVertices[i]);
        Edge* sideEdge2 = new Edge(bottomVertices[next], topVertices[next]);

        // Horizontale Kanten
        Edge* bottomEdge = new Edge(bottomVertices[i], bottomVertices[next]);
        Edge* topEdge = new Edge(topVertices[i], topVertices[next]);

        solid->addEdge(sideEdge1);
        solid->addEdge(sideEdge2);
        solid->addEdge(bottomEdge);
        solid->addEdge(topEdge);

        // Seitenfläche
        vector<Edge*> sideFaceEdges = { sideEdge1, topEdge, sideEdge2, bottomEdge };
        Face* sideFace = new Face(sideFaceEdges);

        // Normale berechnen (nach außen)
        Vector3D v1(bottomVertices[next]->getPosition() - bottomVertices[i]->getPosition());
        Vector3D v2(topVertices[i]->getPosition() - bottomVertices[i]->getPosition());
        Vector3D normal = v1.cross(v2).normalized();
        sideFace->setNormal(normal);

        solid->addFace(sideFace);
    }

    // Bodenfläche (Polygon)
    if (!bothSides || true) { // Immer Boden erstellen
        vector<Edge*> bottomFaceEdges;
        for (int i = 0; i < n; i++) {
            int prev = (i == 0) ? n - 1 : i - 1;
            Edge* edge = new Edge(bottomVertices[prev], bottomVertices[i]);
            solid->addEdge(edge);
            bottomFaceEdges.push_back(edge);
        }
        Face* bottomFace = new Face(bottomFaceEdges);
        bottomFace->setNormal(Vector3D(0, 0, -1));
        solid->addFace(bottomFace);
    }

    // Deckfläche (Polygon)
    vector<Edge*> topFaceEdges;
    for (int i = 0; i < n; i++) {
        int next = (i + 1) % n;
        Edge* edge = new Edge(topVertices[i], topVertices[next]);
        solid->addEdge(edge);
        topFaceEdges.push_back(edge);
    }
    Face* topFace = new Face(topFaceEdges);
    topFace->setNormal(Vector3D(0, 0, 1));
    solid->addFace(topFace);

    // Topologie verlinken
    for (Face* face : solid->getFaces()) {
        for (Edge* edge : face->getEdges()) {
            edge->addFace(face);
        }
    }

    for (Edge* edge : solid->getEdges()) {
        edge->getV1()->addEdge(edge);
        edge->getV2()->addEdge(edge);
    }

    solid->updateBoundingBox();
    return solid;
}

Solid* SketcherToSolid::revolveSketch(const Sketch* sketch, double angleDegrees, int segments) {
    if (!sketch || segments < 3) {
        throw invalid_argument("Invalid parameters");
    }

    vector<Point2D> profile = sketch->getClosedProfile();
    if (profile.size() < 2) {
        throw runtime_error("Sketch profile too small");
    }

    double angleRad = angleDegrees * M_PI / 180.0;
    if (angleRad <= 0 || angleRad > 2 * M_PI) {
        throw invalid_argument("Invalid revolution angle");
    }

    Solid* solid = new Solid();

    // Vertices erstellen
    vector<vector<Vertex*>> vertexGrid;

    for (size_t i = 0; i < profile.size(); i++) {
        vector<Vertex*> segmentVertices;
        double y = profile[i].y;  // Revolutionsachse ist Y-Achse
        double radius = abs(profile[i].x);  // Abstand von der Achse

        for (int j = 0; j <= segments; j++) {
            double theta = angleRad * j / segments;
            double x = radius * cos(theta);
            double z = radius * sin(theta);

            // Für negative x-Werte, Normale umdrehen
            if (profile[i].x < 0) {
                x = -x;
                z = -z;
            }

            Vertex* v = new Vertex(Point3D(x, y, z));
            segmentVertices.push_back(v);
            solid->addVertex(v);
        }
        vertexGrid.push_back(segmentVertices);
    }

    // Faces erstellen (Quads)
    for (size_t i = 0; i < profile.size() - 1; i++) {
        for (int j = 0; j < segments; j++) {
            Vertex* v00 = vertexGrid[i][j];
            Vertex* v01 = vertexGrid[i][j + 1];
            Vertex* v10 = vertexGrid[i + 1][j];
            Vertex* v11 = vertexGrid[i + 1][j + 1];

            // Kanten
            Edge* e1 = new Edge(v00, v01);
            Edge* e2 = new Edge(v01, v11);
            Edge* e3 = new Edge(v11, v10);
            Edge* e4 = new Edge(v10, v00);

            solid->addEdge(e1);
            solid->addEdge(e2);
            solid->addEdge(e3);
            solid->addEdge(e4);

            // Face
            vector<Edge*> faceEdges = { e1, e2, e3, e4 };
            Face* face = new Face(faceEdges);

            // Normale berechnen
            Vector3D v01_v00 = v01->getPosition() - v00->getPosition();
            Vector3D v10_v00 = v10->getPosition() - v00->getPosition();
            Vector3D normal = v01_v00.cross(v10_v00).normalized();
            face->setNormal(normal);

            solid->addFace(face);
        }
    }

    // Endkappen (wenn vollständige Revolution)
    if (abs(angleRad - 2 * M_PI) < 0.001) {
        // Startkappe
        vector<Edge*> startCapEdges;
        for (size_t i = 0; i < profile.size() - 1; i++) {
            Edge* edge = new Edge(vertexGrid[i][0], vertexGrid[i + 1][0]);
            solid->addEdge(edge);
            startCapEdges.push_back(edge);
        }
        Face* startCap = new Face(startCapEdges);
        startCap->setNormal(Vector3D(-1, 0, 0));
        solid->addFace(startCap);

        // Endkappe
        vector<Edge*> endCapEdges;
        for (size_t i = 0; i < profile.size() - 1; i++) {
            Edge* edge = new Edge(vertexGrid[i][segments], vertexGrid[i + 1][segments]);
            solid->addEdge(edge);
            endCapEdges.push_back(edge);
        }
        Face* endCap = new Face(endCapEdges);
        endCap->setNormal(Vector3D(1, 0, 0));
        solid->addFace(endCap);
    }

    // Topologie verlinken
    for (Face* face : solid->getFaces()) {
        for (Edge* edge : face->getEdges()) {
            edge->addFace(face);
        }
    }

    for (Edge* edge : solid->getEdges()) {
        edge->getV1()->addEdge(edge);
        edge->getV2()->addEdge(edge);
    }

    solid->updateBoundingBox();
    return solid;
}

Solid* SketcherToSolid::sweepSketch(const Sketch* sketch, const std::vector<Point3D>& path, bool followRotation) {
    // Fortgeschrittene Funktion für später
    throw runtime_error("Sweep not implemented yet");
}