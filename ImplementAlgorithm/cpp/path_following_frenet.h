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
/* for plotting */
#include "lib/matplotlibcpp.h"

class pathFollowingFrenet 
{
  private:
    /* fixed parameters */
    const double dt = 0.05;                  // time step
    const int    N  = 30;                    // iteration steps
    const double vr = 10.0;                  // longditudinal velocity

    /* state variables */ 
    double s_;
    double e_;
    double theta_e_;
    std::vector<double> s_vec_;
    std::vector<double> e_vec_;
    std::vector<double> theta_e_vec_;

    /* store traj of catesian */
    std::vector<double> x0_vec_;
    std::vector<double> y0_vec_;
    std::vector<double> x_vec_;
    std::vector<double> y_vec_;

    /* spline variables */
    NaturalSpline<QVector2D>* spline_ptr_;   // the spline function, spline_ := s(t)
    std::vector<double>* spline_seq_ptr_;    // spline seq, s_i(t), i from 0 to sz_
    double kappa_s_;                         // curvature w.r.t. s, kappa_s_ := kappa(s)
    int sz_;                                 // size of givin points to the spline

    /* control variables */
    double omega_;                           // calculated control input
    double k_theta_e_;
    double k_e_;

    void fitSpline(std::vector<double>& posX, 
                   std::vector<double>& posY);
    void setControlGains(double k_theta_e, 
                         double k_e);
    double reverseArclength();               // obtain t from s(t) by inversing the function
    double calKappa(double t);
    double calOmega();
    void Frenet2Cartesian(double t);
    void propagate();                        // execute the propagation for one step
    void augmentStateVectors();
    void investigateStates() const;
    void plotStates() const;
    void plotXYCartesian() const;

  public:
    pathFollowingFrenet(double s, 
                        double e, 
                        double theta_e);
    ~pathFollowingFrenet();

    void investigateSpline() const;

};

#endif /* PATH_FOLLOWING_FRENET_H */