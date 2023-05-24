function [ T, rp, ra, eps, b ] = orbitParameters( a, e, i, OM, om, theta, mu)
% Calculates main orbit parameters
%
% PROTOTYPE:
%  [ T, rp, ra, eps, b ] = orbitParameters( a, e, i, OM, om, theta, mu)
%
% INPUTS:
%   a    [1]    Semi-major axis        [km]
%   e    [1]    Eccentricity           [-]
%   i    [1]    Inclination            [rad]
%   OM   [1]    RAAN                   [rad]
%   om   [1]    Argument of pericenter [rad]
%   theta[1]    Real anomaly           [rad]
%   mu   [1]    Planetary constant     [km^3/s^2]
%
% OUTPUTS:
%   T    [1]    Period of orbit        [s]
%   rp   [1]    Radius of pericenter   [km]
%   ra   [1]    Radius of apocenter    [km]
%   eps  [1]    Specific energy        [km^2/s^2]
%   b    [1]    Semi-minor axis        [km]
%                   -For hyperbolic orbits, b = b_bar
%
% CONTRIBUTORS:
%   Matteo D'Ambrosio

T = sqrt( 4 .* pi.^2 .* a.^3 ./ mu ); % [s] - Period of orbit

eps = -mu ./ ( 2 .* a ); % [J/kg] - Specific energy

[ r, v ] = kp2rv( a, e, i, OM, om, theta, mu );

h = cross( r , v );
h = norm(h);

P = h.^2 ./ mu;

rp = P ./ (1 + e .* cos(0) ); % [km] - Radius at pericenter
ra = P ./ (1 + e .* cos(pi) ); % [km] - Radius at apocenter

b = abs(a) .* sqrt( abs( 1 - e.^2 ) );

end