#include "lib/spline_library/splines/natural_spline.h"
#include "lib/spline_library/vector.h"
#include "lib/spline_library/spline.h"
#include <iostream>
#include <QVector2D>

int main(int argc, char** argv)
{   
    // I guess:
    // we need to input P_{i,t} = ( x_i(t), y_i(t) )
    // and curvature can give us d^2 P_{i,t}/dt^2 = ( x_i''(t), y_i''(t) )
    std::vector<QVector2D> splinePoints{
        QVector2D( 0, 0),
        QVector2D( 1, 3),
        QVector2D( 4, 2),
        QVector2D( 6, 4),
        QVector2D( 8, 4),
        QVector2D( 7, 5)
    };

    NaturalSpline<QVector2D> mySpline(splinePoints);
    QVector2D p0 = mySpline.getPosition(0.0f);
    QVector2D p1 = mySpline.getPosition(1.0f);
    QVector2D p2 = mySpline.getPosition(2.0f);
    QVector2D p3 = mySpline.getPosition(3.0f);
    QVector2D p4 = mySpline.getPosition(4.0f);
    QVector2D p5 = mySpline.getPosition(4.05f);
    double arclength = mySpline.arcLength(4, 5);
    double totallength = mySpline.totalLength();
    auto curvatureSpline = mySpline.getCurvature(1.0f);
    

    std::cout << "arclength = " << arclength
              << ", totallength = " << totallength
              << ", curvature[0] = " << curvatureSpline.curvature[0]
              << ", curvature[1] = " << curvatureSpline.curvature[1]
              << ", p0[0] = " << p0[0]
              << ", p0[1] = " << p0[1] 
              << ", p1[0] = " << p1[0]
              << ", p1[1] = " << p1[1]
              << ", p2[0] = " << p2[0]
              << ", p2[1] = " << p2[1]
              << ", p3[0] = " << p3[0]
              << ", p3[1] = " << p3[1]
              << ", p4[0] = " << p4[0]
              << ", p4[1] = " << p4[1]
              << ", p5[0] = " << p5[0]
              << ", p5[1] = " << p5[1] << std::endl;
    
    return 0;
}