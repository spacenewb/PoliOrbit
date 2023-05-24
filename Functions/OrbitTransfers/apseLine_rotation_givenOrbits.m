function [ dV1, dGamma1, theta11, theta21, dV2, dGamma2, theta12, theta22 ] = apseLine_rotation_givenOrbits( delta_om, rp1, ra1, rp2, ra2, mu )
% Solves apse line rotation maneuver given delta_om, initial and final
% orbits as input
%
% PROTOTYPE:
%   [ dV1, dGamma1, theta11, theta21, dV2, dGamma2, theta12, theta22 ] = apseLine_rotation_givenOrbits( delta_om, rp1, ra1, rp2, ra2, mu )
%
% INPUTS:
%   delta_om[1]     Apse line rotation                          [rad]
%   rp1[1]          Orbit 1 pericenter radius                   [km]
%   ra1[1]          Orbit 1 apocenter radius                    [km]
%   rp2[1]          Orbit 2 pericenter radius                   [km]
%   ra2[1]          Orbit 2 apocenter radius                    [km]
%   mup[1]          Planetary constant                          [km^3/s^2]
%
% OUTPUTS:
%   dV1[1]           deltaV required at solution 1              [km/s]
%   dGamma1[1]       Flight path angle of deltaV at solution 1  [rad]
%   theta11[1]       True anomaly of solution 1 on orbit 1      [rad]
%   theta21[1]       True anomaly of solution 1 on orbit 2      [rad]
%   dV2[1]           deltaV required at solution 2              [km/s]
%   dGamma2[1]       Flight path angle of deltaV at solution 2  [rad]
%   theta12[1]       True anomaly of solution 2 on orbit 1      [rad]
%   theta22[1]       True anomaly of solution 2 on orbit 2      [rad]

[ e1, h1 ] = calc_eh( rp1, ra1, mu );
[ e2, h2 ] = calc_eh( rp2, ra2, mu );

A = h2.^2 .* e1 - h1.^2 .* e2;
B = -h1.^2 .* e2 .* sin( delta_om );
C = h1.^2 - h2.^2;

phi = atan( B / A );

% Solution 1
theta11 = phi + acos( C / A * cos( phi ) );
theta21 = theta11 - delta_om;
r_I = ( h1.^2 ./ mu ) ./ ( 1 + e1 .* cos(theta11) );
v_t1 = h1 / r_I;
v_r1 = ( mu / h1 ) .* e1 .* sin(theta11);
v_t2 = h2/r_I;
v_r2 = ( mu / h2 ) .* e2 .* sin(theta21);
dVr = v_r2 - v_r1;
dVt = v_t2 - v_t1;
dV = [ dVr; dVt ];
dV1 = norm( dV );
dGamma1 = deltaGamma( dVr, dVt );

% Solution 2
theta12 = phi - acos( C / A * cos( phi ) );
theta22 = theta12 - delta_om;
r_J = ( h1.^2 ./ mu ) ./ ( 1 + e1 .* cos(theta12) );
v_t1 = h1 / r_J;
v_r1 = ( mu / h1 ) .* e1 .* sin(theta12);
v_t2 = h2/r_J;
v_r2 = ( mu / h2 ) .* e2 .* sin(theta22);
dVr = v_r2 - v_r1;
dVt = v_t2 - v_t1;
dV = [ dVr; dVt ];
dV2 = norm( dV );
dGamma2 = deltaGamma( dVr, dVt );

end

%% Local functions

function [ e, h ] = calc_eh( rp, ra, mu )

e = ( ra - rp ) / ( ra + rp ) ;
h = sqrt( 2 .* mu * ra * rp / ( ra + rp) );

end

function dGamma = deltaGamma( dVr, dVt )

if dVt > 0

    dGamma = atan( dVr / dVt );

else

    dGamma = pi - abs( atan( dVr / dVt ) );

end

end

