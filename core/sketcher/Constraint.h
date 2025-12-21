#pragma once

enum class ConstraintType {
    Horizontal,
    Vertical,
    Distance,
    Angle,
    Tangent,
    Coincident  
};

struct Constraint {
    ConstraintType type;

    // Indizes auf Curves oder Points
    int entityA = -1;
    int entityB = -1;

    double value = 0.0; // Abstand, Winkel, etc.
};
