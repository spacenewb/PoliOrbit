function dy = ode_2bp( ~, y, mu, pert)

% ode_2bp   ODE system for the two body problem ( Keplerian motion)
%
% PROTOTYPE
%   dy = ode_2bp( t, y, mu, pert)
%
% INPUT:
%   t[1]        Time (can be omitted, as the system is autonomous) [T]
%   y[6x1]      State of the body ( rx , ry , rz , vx , vy , vz ) [ L, L/T ]
%   mu[1]       Gravitational parameter of the primary [ L^3/T^2 ]
%   pert[-]     String describing what perturbation to implement ['J2'/'0']
%
% OUTPUT:
%   dy[6x1]     Derivative of the state [ L/T^2, L/T^3 ]
%
% CALLED FUNCTIONS:
%   astroConstants.m
%
% CONTRIBUTORS:
%   Hariharan Venkatesh Vitaladevuni
%
% VERSIONS:
%   2021-10-03: First version
%
%--------------------------------------------------------------------------

% Position and velocity
r = y(1:3);
v = y(4:6);

% Distance from the primary
rnorm = norm(r)

% Calculate perturbation acceleration
switch pert
    case '0'
        acc_pert = [0;0;0];
    case 'J2'
        J2 = astroConstants(9); % Radius of earth value compatibility check
        R_E = 6378.1363; % [km] Check documentation of astroConstants
        J2_const_1 = 3/2*J2*mu*R_E^2/rnorm^4;
        J2_const_2 = r/rnorm.*[(5*r(3)^2/rnorm^2 - 1); (5*r(3)^2/rnorm^2 - 1); (5*r(3)^2/rnorm^2 - 3)];
        acc_pert = J2_const_1*J2_const_2
end

% Set the derivatives of the state
dy = [ v; (-mu/rnorm^3)*r + acc_pert ];

end