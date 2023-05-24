function [ r, v ] = kp2rv( a, e, i, OM, om, theta, mu )
% kp2rv Calculates r,v from Keplerian parameters
%
% PROTOTYPE: 
%   [ r, v ] = kp2rv( a, e_norm, i, OM, om, theta, mu )
%
% INPUT:     SIZE:        DESCRIPTION:        UNITS:
%   a          [1x1]        Semi-axis major     [km]
%   e          [1x1]        Eccentricity        [-]
%   i          [1x1]        Inclination         [rad]
%   OM         [1x1]        RAAN                [rad]
%   om         [1x1]        Pericenter anomaly  [rad]
%   theta      [1x1]        True anomaly        [rad]
%   mu         [1x1]        Planetary constant  [km^3/s^2]
%
% OUTPUT:     SIZE:        DESCRIPTION:        UNITS:
%   r          [3x1]        Position vector     [km]
%   v          [3x1]        Velocity vector     [km/s]
%

%% Calculate r, v in perifocal frame

unit_r = [cos(theta); sin(theta); 0];
unit_theta = [-sin(theta); cos(theta); 0];
P = a.*(1-e.^2);

r_pf = a .* (1 - e.^2) ./ (1 + e .* cos(theta)) * unit_r; % r in perifocal plane
v_r = sqrt( mu ./ P ) .* e .* sin(theta) * unit_r;
v_theta = sqrt( mu ./ P ) .* ( 1 + e .* cos(theta) ) * unit_theta;

v_pf = v_r + v_theta;

%% Generate rotation matrix

R_ECEI2PF = ECEI2pf_rotate( i, OM, om );
r = R_ECEI2PF'*r_pf;
v = R_ECEI2PF'*v_pf;
