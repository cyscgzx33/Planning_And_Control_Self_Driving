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
    delete spline_seq_ptr_;
}

void pathFollowingFrenet::setControlGains(double k_theta_e, double k_e)
{
    k_theta_e_ =  k_theta_e;
    k_e_ = k_e;
}

double pathFollowingFrenet::calKappa()
{
    // given s := s(t), need to solve kappa(s(t))
    // ==> solve t first, and obtain kappa(s(t)) correspondingly

    /* using binary search to resolve value of t from s(t) */
    // firstly find the lowerbound ( TC: O(logN) )
    auto lb = std::lower_bound(spline_seq_ptr_->begin(), spline_seq_ptr_->end(), s_);
    assert ( lb == spline_seq_ptr_->end() ); // all elements in spline_seq_ptr_ is larger than s_
    assert ( lb == spline_seq_ptr_->begin() );
    int lb_idx = lb - spline_seq_ptr_->begin();

    // TODO: find a good way to handle the boundary probelm
    // if (lb_idx == 0 || lb_idx == sz_)
    //     return double(lb_idx);
    
    // start binary search
    double tolerance = 0.000001; // 10^-6
    double start = lb_idx - 1, end = lb_idx;
    while ( end - start > tolerance)
    {
        double mid = start + (end - start) / 2;
        if (spline_ptr_->arcLength(0, mid) < s_)
            start = mid;
        else
            end = mid; 
    }

    /* calculate kappa using curvature formula */
    double t = start;

}

void pathFollowingFrenet::calOmega()
{
    // execute a simple feedback control law to calculate control input
    omega_ = vr * kappa_s_ * cos(theta_e_) / (1 - kappa_s_ * e_) - k_theta_e_* vr * theta_e_ 
             - k_e_ * vr * sin(theta_e_) / theta_e_ * e_;
}

void pathFollowingFrenet::fitSpline(std::vector<double>& posX, std::vector<double>& posY)
{   
    sz_ = posX.size();

    std::vector<QVector2D> spline_pts;
    for (int i = 0; i < sz_; i++)
        spline_pts.push_back( QVector2D( posX[i], posY[i] ) );

    // Ethan: not optimal to use pointer to init a NaturalSpline object here
    // TODO: figure out if another without using pointer would work
    NaturalSpline<QVector2D>* spline_ptr = new NaturalSpline<QVector2D>(spline_pts);
    spline_ptr_ = spline_ptr;

    // TODO: figure out if another without using pointer would work
    std::vector<double>* spline_seq_ptr = new std::vector<double>(sz_);
    for ( int i = 0; i < sz_; i++ )
        (*spline_seq_ptr)[i] = spline_ptr->arcLength(0, i);
    spline_ptr_ = spline_ptr;
}

void pathFollowingFrenet::propagate()
{   
    // update the curvature kappa_s_
    calKappa();
    
    // obtain the control input omega_
    calOmega();

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