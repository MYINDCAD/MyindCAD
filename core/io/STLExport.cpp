// STLExport.cpp
#include "STLExport.h"
#include "solid.h"
#include "Face.h"
#include "Vertex.h"
#include "Point3D.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>

using namespace std;

namespace {
    // Hilfsfunktion zum Schreiben von float im IEEE 754 Format
    void writeFloat(ostream& os, float value) {
        os.write(reinterpret_cast<const char*>(&value), sizeof(float));
    }

    // Hilfsfunktion zum Schreiben von Vector3D
    void writeVector(ostream& os, const Vector3D& vec) {
        writeFloat(os, static_cast<float>(vec.x));
        writeFloat(os, static_cast<float>(vec.y));
        writeFloat(os, static_cast<float>(vec.z));
    }

    // Dreieck triangulieren (für nicht-dreieckige Faces)
    vector<Vector3D> triangulateFace(const Face* face) {
        vector<Vector3D> triangles;

        // Einfache Triangulation: Fan von erstem Vertex
        if (face->getVertices().size() >= 3) {
            const vector<Vertex*>& vertices = face->getVertices();
            Vector3D v0 = vertices[0]->getPosition();

            for (size_t i = 1; i < vertices.size() - 1; i++) {
                triangles.push_back(v0);
                triangles.push_back(vertices[i]->getPosition());
                triangles.push_back(vertices[i + 1]->getPosition());
            }
        }

        return triangles;
    }
}

bool STLExport::exportToSTL(const Solid* solid, const string& filename, bool binary) {
    if (!solid || solid->getFaces().empty()) {
        cerr << "Error: Solid is empty or invalid" << endl;
        return false;
    }

    if (binary) {
        return exportToSTLBinary(solid, filename);
    }
    else {
        return exportToSTLAscii(solid, filename);
    }
}

bool STLExport::exportToSTLBinary(const Solid* solid, const string& filename) {
    ofstream file(filename, ios::binary);
    if (!file) {
        cerr << "Error: Cannot open file " << filename << endl;
        return false;
    }

    // Header (80 Bytes)
    char header[80];
    strncpy(header, "MyIndustrialCAD STL Export", 79);
    header[79] = '\0';
    file.write(header, 80);

    // Anzahl der Dreiecke zählen
    uint32_t triangleCount = 0;
    for (const Face* face : solid->getFaces()) {
        int verticesCount = face->getVertices().size();
        if (verticesCount >= 3) {
            triangleCount += (verticesCount - 2); // n-Eck hat n-2 Dreiecke
        }
    }

    file.write(reinterpret_cast<const char*>(&triangleCount), sizeof(uint32_t));

    // Dreiecke schreiben
    for (const Face* face : solid->getFaces()) {
        const vector<Vertex*>& vertices = face->getVertices();
        if (vertices.size() < 3) continue;

        Vector3D normal = face->getNormal();

        // Einfache Triangulation: Triangle Fan
        for (size_t i = 1; i < vertices.size() - 1; i++) {
            // Normale
            writeVector(file, normal);

            // Vertex 0
            writeVector(file, vertices[0]->getPosition());

            // Vertex i
            writeVector(file, vertices[i]->getPosition());

            // Vertex i+1
            writeVector(file, vertices[i + 1]->getPosition());

            // Attribute byte count (0)
            uint16_t attribute = 0;
            file.write(reinterpret_cast<const char*>(&attribute), sizeof(uint16_t));
        }
    }

    file.close();
    cout << "Exported " << triangleCount << " triangles to " << filename << endl;
    return true;
}

bool STLExport::exportToSTLAscii(const Solid* solid, const string& filename) {
    ofstream file(filename);
    if (!file) {
        cerr << "Error: Cannot open file " << filename << endl;
        return false;
    }

    file << "solid " << solid->getName() << endl;

    int triangleCount = 0;
    for (const Face* face : solid->getFaces()) {
        const vector<Vertex*>& vertices = face->getVertices();
        if (vertices.size() < 3) continue;

        Vector3D normal = face->getNormal();

        for (size_t i = 1; i < vertices.size() - 1; i++) {
            file << "  facet normal "
                << normal.x << " " << normal.y << " " << normal.z << endl;
            file << "    outer loop" << endl;

            file << "      vertex "
                << vertices[0]->getPosition().x << " "
                << vertices[0]->getPosition().y << " "
                << vertices[0]->getPosition().z << endl;

            file << "      vertex "
                << vertices[i]->getPosition().x << " "
                << vertices[i]->getPosition().y << " "
                << vertices[i]->getPosition().z << endl;

            file << "      vertex "
                << vertices[i + 1]->getPosition().x << " "
                << vertices[i + 1]->getPosition().y << " "
                << vertices[i + 1]->getPosition().z << endl;

            file << "    endloop" << endl;
            file << "  endfacet" << endl;

            triangleCount++;
        }
    }

    file << "endsolid " << solid->getName() << endl;
    file.close();

    cout << "Exported " << triangleCount << " triangles to " << filename << " (ASCII)" << endl;
    return true;
}

bool STLExport::exportToOBJ(const Solid* solid, const string& filename) {
    ofstream file(filename);
    if (!file) {
        cerr << "Error: Cannot open file " << filename << endl;
        return false;
    }

    file << "# MyIndustrialCAD OBJ Export" << endl;
    file << "# Vertices: " << solid->getVertices().size() << endl;
    file << "# Faces: " << solid->getFaces().size() << endl;
    file << endl;

    // Vertices
    int vertexIndex = 1;
    for (const Vertex* vertex : solid->getVertices()) {
        const Point3D& p = vertex->getPosition();
        file << "v " << p.x << " " << p.y << " " << p.z << endl;
    }

    file << endl;

    // Face Normals (optional)
    file << "# Face normals" << endl;
    for (const Face* face : solid->getFaces()) {
        const Vector3D& n = face->getNormal();
        file << "vn " << n.x << " " << n.y << " " << n.z << endl;
    }

    file << endl << "# Faces" << endl;

    // Faces mit Vertex-Indizes
    int faceIndex = 0;
    for (const Face* face : solid->getFaces()) {
        const vector<Vertex*>& vertices = face->getVertices();
        if (vertices.empty()) continue;

        file << "f ";
        for (size_t i = 0; i < vertices.size(); i++) {
            // Finde Vertex Index (1-basiert in OBJ)
            int vIndex = 1;
            const auto& allVertices = solid->getVertices();
            for (size_t j = 0; j < allVertices.size(); j++) {
                if (allVertices[j] == vertices[i]) {
                    vIndex = static_cast<int>(j) + 1;
                    break;
                }
            }

            // Vertex/Texture/Normal Index (wir haben keine Texturen)
            file << vIndex << "//" << (faceIndex + 1);
            if (i < vertices.size() - 1) file << " ";
        }
        file << endl;

        faceIndex++;
    }

    file.close();
    cout << "Exported to OBJ: " << filename << endl;
    return true;
}