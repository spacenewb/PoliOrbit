function [eps,h_vec,e,cos_f,v_r,v_t,e_dot_h] = orb_mot_const(Y, mu)

%orb_mot_const Calculate evolution of orbital constants 
%               for the two body problem ( Keplerian motion). Function is
%               vectorised, avoid for-loops for performance.
%
% PROTOTYPE:
%   [eps,h,e,cos_f,v_r,v_t,e_dot_h] = orb_mot_const(Y,mu)
%
% INPUT:
%   Y[6x1]      State of the body ( rx , ry , rz , vx , vy , vz ) [ L, L/T ]
%   mu[1]       Gravitational parameter of the primary [ L^3/T^2 ]
%
% OUTPUT:
%   T[1]        Orbital period [ T ]
%   eps[1]      Specific orbital energy [ L^2/T^2 ]
%   h[3x1]      Specific angular momentum vector []
%   e[3x1]      Eccentricity vector [-]
%   cos_f[1]    Cosine of true anomaly [-]
%   v_r[1]      Radial velocity [ L/T ]
%   v_t[1]      Transversal velocity [ L/T ]
%
% CALLED FUNCTIONS:
%   [N.A]
%
% CONTRIBUTORS:
%   Hariharan Venkatesh Vitaladevuni
%
% VERSIONS:
%   2021-11-13: First version
%
%--------------------------------------------------------------------------

% Restructuring Input
r_vec = Y(1:3,:);
v_vec = Y(4:6,:);

% Calculation
%a = 1./(2./vecnorm(r) - dot(v,v,1)./mu); % Semi major axis [km]
%T = 2.*pi.*(a.^3./mu).^0.5; % Time period of Orbit [s]

eps = vecnorm(v_vec).^2./2 - mu./vecnorm(r_vec); % Specific energy 
h_vec = cross(r_vec,v_vec,1); % Angular momentum
e = cross(v_vec./mu, h_vec, 1) - r_vec./vecnorm(r_vec); % Eccentricity
cos_f = dot(r_vec,e,1)./(vecnorm(r_vec).*vecnorm(e)); % Cosine of True Anomaly

v_r = dot(v_vec,r_vec./vecnorm(r_vec),1); % Radial velocity [km/s]
v_t = (vecnorm(v_vec).^2 - v_r.^2).^0.5; % Transversal velocity [km/s]

e_dot_h = dot(e, h_vec, 1); % Perpendicularity of e and h

end
