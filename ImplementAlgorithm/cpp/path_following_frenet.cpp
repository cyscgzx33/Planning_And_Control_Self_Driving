#include "path_following_frenet.h"

pathFollowingFrenet::pathFollowingFrenet(double s, double e, double theta_e) : 
                     s_(s), e_(e), theta_e_(theta_e), k_theta_e_(1.0), k_e_(1.0) 
{
    // initialize the spline
    std::vector<double> posX = {1.0, 2.0, 3.0};
    std::vector<double> posY = {2,0, 3.0, 6.0};
    fitSpline(posX, posY);
    std::cout << "successfully inited a pathFollowingFrenet object" << std::endl;
}

pathFollowingFrenet::~pathFollowingFrenet() 
{
    delete spline_ptr_;
}

void pathFollowingFrenet::setControlGains(double k_theta_e, double k_e)
{
    k_theta_e_ =  k_theta_e;
    k_e_ = k_e;
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
    std::vector<QVector2D> spline_pts;
    for (int i = 0; i < posX.size(); i++)
        spline_pts.push_back( QVector2D( posX[i], posY[i] ) );

    // Ethan: not optimal to use pointer to init a NaturalSpline object here
    // TODO: figure out if another without using pointer would work
    NaturalSpline<QVector2D>* spline_ptr = new NaturalSpline<QVector2D>(spline_pts);
    spline_ptr_ = spline_ptr;
}

void pathFollowingFrenet::propagate()
{
    // execute the kinematics equation to update the state variables
    s_        =  s_ + vr * cos(theta_e_) / (1 - kappa_s_ * e_) * dt;
    e_        =  e_ + vr * sin(theta_e_) * dt;
    theta_e_  =  theta_e_ + ( omega_ - vr * kappa_s_ * cos(theta_e_) / (1 - kappa_s_ * e_) ) * dt;
}

void pathFollowingFrenet::investigateSpline() const
{
    std::cout << "arclength = " << spline_ptr_->totalLength() << std::endl;
}

int main(int argc, char** argv)
{
    pathFollowingFrenet scenario_0(0.1, 1.0, 0.05);
    scenario_0.investigateSpline();

    return 0;
}