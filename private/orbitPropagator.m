function [tspan,Y] = orbitPropagator(y0, mu, pert, n_orbs, n_eval)

% orbitPropagator   ODE system for the two body problem ( Keplerian motion)
%
% PROTOTYPE
%   dy = ode_2bp( t, y, mu, pert)
%
% INPUT:
%   y0[6x1]     Initial state of the body ( rx , ry , rz , vx , vy , vz ) [ L, L/T ]
%   mu[1]       Gravitational parameter of the primary [ L^3/T^2 ]
%   pert[-]     String describing what perturbation to implement ['J2'/'0']
%   n_orbs[1]   Nmber of orbits to propagate [ - ]
%   n_eval[1]   Number of points to evaluate the propagation per orbit [ - ]
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

if isempty(pert)
    pert = [];
end

%% Restructuring Input
r0 = y0(1:3)';
v0 = y0(4:6)';

%% Set time span
a = 1 /( 2 /norm(r0) - dot(v0,v0)/ mu ); % Semi major axis [km]
T = 2*pi*sqrt( a^3/mu ); % Orbital period [1/s]
tspan = linspace ( 0, n_orbs*T, n_eval*n_orbs);

%% Orbit Propagation
% Set options for the ODE solver
options = odeset ( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[ ~, Y ] = ode113( @(t,y) ode_2bp_perturbed(t,y,mu,pert), tspan , y0, options );
Y = Y';

end