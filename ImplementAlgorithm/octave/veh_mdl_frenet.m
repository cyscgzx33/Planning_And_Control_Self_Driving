function [s, e, theta_e] = veh_mdl_frenet(s, e, theta_e)

% A fixed time step dt = 0.05 s
dt = 0.05;

% A fixed longitudinal velocity vr as 10 m/s 
vr = 10;

% obtain the kappa_s according to s(t)


% Propagate the formula
s        =  s + vr * cos(theta_e) / (1 - kappa_s * e) * dt;
e        =  e + vr * sin(theta_e) * dt;
theta_e  =  theta_e + ( omega - vr * kappa_s * cos(theta_e) / (1 - kappa_s * e) ) * dt;

end