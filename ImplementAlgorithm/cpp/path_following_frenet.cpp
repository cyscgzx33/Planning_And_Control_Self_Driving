#include "path_following_frenet.h"

namespace plt = matplotlibcpp;

pathFollowingFrenet::pathFollowingFrenet(double s, double e, double theta_e) : 
                     s_(s), e_(e), theta_e_(theta_e), k_theta_e_(1.0), k_e_(1.0) 
{
    /* initialize the spline */
    std::vector<double> posX = {0.0, 1.0, 3.0, 4.0, 5.0, 5.5, 6.5, 7.7, 9.0,  10.5, 12,   14,   16,   10};
    std::vector<double> posY = {0,0, 2.0, 4.0, 4.5, 6.5, 8.1, 9.4, 7.5, 10.0, 11.5, 12.5, 14.6, 16.9, 24};
    fitSpline(posX, posY);

    /* start path following */
    std::cout << "start to execute the simulation of propagation" << std::endl;
    std::cout << "The initial vehicle states are:" << std::endl;
    investigateStates();
    augmentStateVectors();
    for (int i = 0; i < N; i++)
    {
        std::cout << "timestep @ i = " << i << ": ";
        propagate();
        investigateStates();
        augmentStateVectors();
    }

    /* plot the resolved state vectors */
    plotStates();

    /* plot paths in Cartesian coordinate */
    plotXYCartesian();
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

double pathFollowingFrenet::reverseArclength()
{
    /** Goal:    given s := s(t), need to solve kappa(s(t))
     *  Method:  solve t first, and obtain kappa(s(t)) correspondingly
     * */

    /* using binary search to resolve value of t from s(t) */
    // firstly find the lowerbound ( TC: O(logN) )
    auto lb = std::lower_bound(spline_seq_ptr_->begin(), spline_seq_ptr_->end(), s_);
    assert ( lb != spline_seq_ptr_->end() ); // all elements in spline_seq_ptr_ is larger than s_

    int lb_idx = lb - spline_seq_ptr_->begin();

    // TODO: find a good way to handle the boundary issue
    if (lb_idx == 0)
        return double (lb_idx);

    // Binary search ( TC: O(logN) )
    double tolerance = 0.000001; // 10^-6
    double start = lb_idx - 1, end = lb_idx;
    while (end - start > tolerance)
    {
        double mid = start + (end - start) / 2;
        if (spline_ptr_->arcLength(0, mid) < s_)
            start = mid;
        else
            end = mid; 
    }

    // start and end are very close to each other
    return start; 
}

double pathFollowingFrenet::calKappa(double t)
{
    /* calculate kappa using curvature formula */
    auto d_spline = spline_ptr_->getCurvature(t); // derivatives of spline @ t, including 1st and 2nd derivatives

    // denote s(t) = ( x(t), y(t) )
    // dsdt(t) = ( dxdt(t), dydt(t) )
    // d2sdt2(t) = ( d2xdt2(t), d2ydt2(t) )
    // according to [url=https://www.math24.net/curvature-radius/]:
    // kappa = abs(dxdt * d2ydt2 - dydt * d2xdt2) / pow(dxdt * dxdt + dydt * dydt, 1.5)
    double dxdt = d_spline.tangent[0];
    double dydt = d_spline.tangent[1];
    double d2xdt2 = d_spline.curvature[0];
    double d2ydt2 = d_spline.curvature[1];
    double kappa_s = (dxdt * d2ydt2 - dydt * d2xdt2) / pow(dxdt * dxdt + dydt * dydt, 1.5);
    
    return kappa_s;
}

double pathFollowingFrenet::calOmega()
{
    // execute a simple feedback control law to calculate control input
    double omega = vr * kappa_s_ * cos(theta_e_) / (1 - kappa_s_ * e_) - k_theta_e_* vr * theta_e_ 
                 - k_e_ * vr * sin(theta_e_) / theta_e_ * e_;

    return omega;
}

void pathFollowingFrenet::Frenet2Cartesian(double t)
{
    /* Step I: obtain the cartesian coordinates on the nominal path */
    auto d_spline = spline_ptr_->getCurvature(t); // derivatives of spline @ t, including 1st and 2nd derivatives
    double x_cart     =  d_spline.position[0];
    double y_cart     =  d_spline.position[1];
    double dxdt_cart  =  d_spline.tangent[0];
    double dydt_cart  =  d_spline.tangent[1];
    // According to [https://math.libretexts.org/Bookshelves/Calculus/Supplemental_Modules_(Calculus)/Vector_Calculus/2%3A_Vector-Valued_Functions_and_Motion_in_Space/2.3%3A_Curvature_and_Normal_Vectors_of_a_Curve]
    // Tangent vector T = dr/dt, where r = x(t)i + y(t)j= ( x(t), y(t) )
    // Therefor, T = (dxdt, dydt)
    // Thus, theta(t) = atan2(dydt, dxdt)
    double theta = atan2(dydt_cart, dxdt_cart);

    /* Step II: obtain the cartesian coordinates of actual position via error equations */
    double x_cart_actual = x_cart + e_ * cos(M_PI + theta + theta_e_);
    double y_cart_actual = y_cart + e_ * sin(M_PI + theta + theta_e_);

    /* Step III: store the results */
    x_vec_.push_back(x_cart_actual);
    y_vec_.push_back(y_cart_actual);
    x0_vec_.push_back(x_cart);
    y0_vec_.push_back(y_cart);
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
    spline_seq_ptr_ = spline_seq_ptr;

    std::cout << "successfully inited a pathFollowingFrenet object" << std::endl;
}

void pathFollowingFrenet::propagate()
{   
    // obtain t by reversing the arclength s := s(t)
    double t = reverseArclength();

    // convert current Frenet pose to Cartesian pose
    Frenet2Cartesian(t);

    // update the curvature kappa_s_
    kappa_s_ = calKappa(t);
    
    // update the control input omega_
    omega_ = calOmega();

    // execute the kinematics equation to update the state variables
    s_        =  s_ + vr * cos(theta_e_) / (1 - kappa_s_ * e_) * dt;
    e_        =  e_ + vr * sin(theta_e_) * dt;
    theta_e_  =  theta_e_ + ( omega_ - vr * kappa_s_ * cos(theta_e_) / (1 - kappa_s_ * e_) ) * dt;
}

void pathFollowingFrenet::augmentStateVectors()
{
    s_vec_.push_back(s_);
    e_vec_.push_back(e_);
    theta_e_vec_.push_back(theta_e_);
}

void pathFollowingFrenet::investigateSpline() const
{
    std::cout << "arclength = " << spline_ptr_->totalLength() << std::endl;
}

void pathFollowingFrenet::investigateStates() const
{
    std::cout << "[s, e, theta_e] = " 
              << s_       << ", "
              << e_       << ", "
              << theta_e_ << "]\n";
}

void pathFollowingFrenet::plotStates() const
{
    /* plot results */
    // It's useful for debugging!
    plt::subplot(3, 1, 1);
    plt::title("arclength: s [m]");
    plt::plot(s_vec_);
    plt::subplot(3, 1, 2);
    plt::title("lateral error: e [m]");
    plt::plot(e_vec_);
    plt::subplot(3, 1, 3);
    plt::title("heading error: theta [rad]");
    plt::plot(theta_e_vec_);

    plt::show();
}

void pathFollowingFrenet::plotXYCartesian() const
{
    /* plot both reference path and actual path in Cartesian coordinates */
    plt::title("Resulst of Path Following Control: (X, Y) [m]");
    plt::named_plot("Actual Path", x_vec_, y_vec_, "b");
    plt::named_plot("Reference Path", x0_vec_, y0_vec_, "r--");
    plt::legend();
    plt::save("./simple_path_following.png");

    plt::show();
}

int main(int argc, char** argv)
{
    pathFollowingFrenet scenario_0(0.1, 1.0, 0.05);
    scenario_0.investigateSpline();
    std::cout << "preparing to destroy object scenario_0\n";

    return 0;
}