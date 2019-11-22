#include <path_following_frenet.h>


pathFollowingFrenet::pathFollowingFrenet(double s, double e, double theta_e) : 
                     s_(s), e_(e), theta_e_(theta_e), k_theta_e_(1.0), k_e_(1.0) {}

pathFollowingFrenet::~pathFollowingFrenet() {}


void pathFollowingFrenet::setControlGains(double k_theta_e, double k_e)
{
    k_theta_e_ =  k_theta_e;
    k_e_ = k_e;

    // initialize the spline
    std::vector<double> posX = {};
    std::vector<double> posY = {};
    fitSpline(posX, posY);
}

void pathFollowingFrenet::calControlInput()
{
    // execute a simple feedback control law to calculate control input
    kappa_s_ = 0.0;
    omega_ = vr * kappa_s_ * cos(theta_e_) / (1 - kappa_s_ * e_) - k_theta_e_* vr * theta_e_ 
             - k_e_ * vr * sin(theta_e_) / theta_e_ * e_;
    
}

void pathFollowingFrenet::fitSpline(std::vector<double>& posX, std::vector<double>& posY)
{
    spline_.set_points(posX, posY);
}

void pathFollowingFrenet::propagate()
{
    // execute the kinematics equation to update the state variables
    s_        =  s_ + vr * cos(theta_e_) / (1 - kappa_s_ * e_) * dt;
    e_        =  e_ + vr * sin(theta_e_) * dt;
    theta_e_  =  theta_e_ + ( omega_ - vr * kappa_s_ * cos(theta_e_) / (1 - kappa_s_ * e_) ) * dt;
}