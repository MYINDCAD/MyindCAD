// STLExport.h
#ifndef STL_EXPORT_H
#define STL_EXPORT_H

#include <string>

class Solid;

class STLExport {
public:
    // Exportiert Solid als STL-Datei
    static bool exportToSTL(const Solid* solid, const std::string& filename, bool binary = true);

    // Exportiert Solid als OBJ-Datei
    static bool exportToOBJ(const Solid* solid, const std::string& filename);

private:
    // Private Implementierungen
    static bool exportToSTLBinary(const Solid* solid, const std::string& filename);
    static bool exportToSTLAscii(const Solid* solid, const std::string& filename);

    STLExport() = delete; // Statische Klasse
};

#endif // STL_EXPORT_H