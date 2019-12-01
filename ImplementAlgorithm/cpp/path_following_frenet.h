/**
 * Author: Yusen Chen
 * 
 * Version
 *  - V0.1, 11/25/2019
 * 
 * Note
 *  - V0.0 (11/19/2019): Implemented a path following algorithm, lateral control only
 *  - V0.1 (11/25/2019): Adding a longitudinal control module
 * */

#ifndef PATH_FOLLOWING_FRENET_H
#define PATH_FOLLOWING_FRENET_H

#include <vector>
#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
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
    const double vr = 10.0;                  // longditudinal velocity (lateral control only)
    const double Kp = 10.0;                  // PID controller params: proporsional
    const double Ki = 0.05;                  // PID controller params: integral
    const double Kd = 2;                     // PID controller params: derivative

    /* state variables */ 
    double s_;
    double e_;
    double theta_e_;
    double s_e_;
    double mu_;                              // stabilizer input for s_e_ / integrator of s_e_
    std::vector<double> s_vec_;
    std::vector<double> e_vec_;
    std::vector<double> theta_e_vec_;

    /* PID controller variables */
    double I_acc_ = 0.0;                     // accumulator variable for integral calculation
    double prev_error_ = 0.0;                // previous term for derivative calculation

    /* store traj of catesian */
    std::vector<double> x0_vec_;
    std::vector<double> y0_vec_;
    std::vector<double> x_vec_;
    std::vector<double> y_vec_;

    /* map information */
    std::vector<std::vector<double>> wps_;   // waypoints read from map csv

    /* spline variables */
    NaturalSpline<QVector2D>* spline_ptr_;   // the spline function, spline_ := s(t)
    NaturalSpline<QVector2D>* v_prf_ptr_;    // reference speed profile spline
    std::vector<double>* spline_seq_ptr_;    // spline seq, s_i(t), i from 0 to sz_
    double kappa_s_;                         // curvature w.r.t. s, kappa_s_ := kappa(s)
    int sz_;                                 // size of givin points to the spline

    /* control variables */
    double omega_;                           // calculated control input
    double k_theta_e_;
    double k_e_;

    void fitSpline(std::vector<double>& posX, 
                   std::vector<double>& posY);
    void fitVelProfile(std::vector<double>& idx,
                       std::vector<double>& vel);
    void setControlGains(double k_theta_e, 
                         double k_e);
    void readRoadmapFromCSV();               // read map info from cvs and store in wps_
    double reverseArclength();               // obtain t from s(t) by inversing the function
    double calKappa(double t);
    double calOmega(double vr_t);            // calculate omega_ via current speed vr_t 
    double calSpeedInput(double t);          // calculate vr_t by calling a PID controller
    void Frenet2Cartesian(double t);
    void propagate();                        // execute the propagation for one step
    void augmentStateVectors();
    void investigateStates() const;
    void plotStates() const;
    void plotXYCartesian() const;

    /* PID controller helper functions */
    double getProp(double error);
    double getIntegral(double error);
    double getDeriv(double error);

  public:
    pathFollowingFrenet(double s, 
                        double e, 
                        double theta_e,
                        double s_e);
    ~pathFollowingFrenet();

    void investigateSpline() const;

};

#endif /* PATH_FOLLOWING_FRENET_H */