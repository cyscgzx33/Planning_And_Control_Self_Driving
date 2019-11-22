#ifndef PATH_FOLLOWING_FRENET_H
#define PATH_FOLLOWING_FRENET_H

#include <vector>
#include <iostream>
#include <cmath>

class pathFollowingFrenet 
{
  private:
    
    /* fixed parameters */
    const double dt = 0.05; // time step
    const double vr = 10.0; // longditudinal velocity

    /* state variables */ 
    double s_;
    double e_;
    double theta_e_;

    /* assistant variables */
    double omega_;
    double kappa_s_;

    /* feedback control gains */
    double k_theta_e_;
    double k_e_;

  public:
    pathFollowingFrenet(double s, double e, double theta_e);
    ~pathFollowingFrenet();

    void fitSpline();
    void setControlGains(double k_theta_e, double k_e);
    void calControlInput();
    void propagate(); // execute the propagation for one step
};

#endif /* PATH_FOLLOWING_FRENET_H */