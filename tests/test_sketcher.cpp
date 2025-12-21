#include "../core/sketcher/Sketch.h"
#include <iostream>

int main() {
    Sketch s;

    // Punkte
    s.addPoint(0, 0);
    s.addPoint(5, 5);
    s.printPoints();

    // Linien
    s.addLine(0, 0, 10, 0);
    s.addLine(0, 0, 0, 5);

    // Kreis & Bogen
    s.addCircle(5, 5, 3);
    s.addArc(5, 5, 3, 0, 90);

    // Rechteck
    s.addRectangle(1, 1, 4, 4);

    // Spline
    s.addSpline({ {0,0},{2,3},{4,1} });

    s.printCurves();

    // Constraints
    s.addHorizontalConstraint(0);
    s.addVerticalConstraint(1);
    s.addDistanceConstraint(0, 1, 5.0);
    s.addAngleConstraint(0, 1, 90.0);

    s.printConstraints();

    return 0;
}
