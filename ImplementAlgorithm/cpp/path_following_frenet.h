#ifndef PATH_FOLLOWING_FRENET_H
#define PATH_FOLLOWING_FRENET_H

#include <vector>
#include <iostream>
#include <cmath>
/* for fitting spline */
#include "lib/spline_library/splines/natural_spline.h"
#include "lib/spline_library/vector.h"
#include "lib/spline_library/spline.h"
#include <QVector2D>

class pathFollowingFrenet 
{
  private:
    
    /* fixed parameters */
    const double dt = 0.05;                  // time step
    const double vr = 10.0;                  // longditudinal velocity

    /* state variables */ 
    double s_;
    double e_;
    double theta_e_;

    /* spline variables */
    NaturalSpline<QVector2D>* spline_ptr_;   // the spline function, spline_ := s(t)
    std::vector<double>* spline_seq_ptr_;    // spline seq, s_i(t), i from 0 to sz_
    double kappa_s_;                         // curvature w.r.t. s, kappa_s_ := kappa(s)
    int sz_;                                 // size of givin points to the spline

    /* control variables */
    double omega_;                           // calculated control input
    double k_theta_e_;
    double k_e_;

  public:
    pathFollowingFrenet(double s, 
                        double e, 
                        double theta_e);
    ~pathFollowingFrenet();

    void fitSpline(std::vector<double>& posX, 
                   std::vector<double>& posY);
    void setControlGains(double k_theta_e, 
                         double k_e);
    double calKappa();
    void calOmega();
    void propagate();                        // execute the propagation for one step
    void investigateSpline() const;
};

#endif /* PATH_FOLLOWING_FRENET_H */