function [ delta_om, e2, h2, theta2 ] = apseLine_rotation_givendeltaV( theta1, deltaV_LVLH, rp1, ra1, r_m, mu )
% Solves apse line rotation maneuver given theta1, deltaV_LVLH as inputs
%
% PROTOTYPE:
%   [ delta_om, e2, h2 ] = apseLine_rotation_givendeltaV( theta1, deltaV_LVLH, rp1, ra1, r_m, mu )
%
% INPUTS:
%   theta1[1]       True anomaly of maneuver on orbit 1                     [rad]
%   deltaV_LVLH[3]  deltaVelocity vector for maneuver given in LVLH frame   [km/s]
%                            LVLH - Radial and transversal directions
%                            deltaV_LVLH = [ dV_r; dV_theta; 0 ];
%   rp1[1]          Radius of perigee of orbit 1                            [km]
%   ra1[1]          Radius of apogee of orbit 1                             [km]
%   r_m[1]          Radius of maneuvering point                             [km]
%   mu[1]           Planetary constant                                      [km^3/s^2]
%
% OUTPUTS:
%   delta_om[1]      Periapsis rotation                         [rad]
%   e2[1]            Eccentricity of second orbit               [-]
%   h2[1]            Angular momentum magnitude of second orbit [km^2/s]
%   theta2[1]        True anomaly on orbit 2                    [rad]
dVr = deltaV_LVLH(1);
dVt = deltaV_LVLH(2);

[ e1, h1 ] = calc_eh( rp1, ra1, mu );
[ v_r1, v_t1 ] = vel_LVLH( e1, h1, theta1, mu );

theta2 = ( v_t1.^2 .* r_m ./ mu ) .* ( v_t1 + dVt ) .* (v_r1 + dVr) ./ ( (v_t1 + dVt).^2 .* e1 .* cos(theta1) + dVt .* ( dVt + 2 .* v_t1 ) ); % True anomaly on second orbit

h2 = h1 + r_m * dVt;
e2 = 1 ./ sin(theta2) .* h2 ./ mu .* ( v_r1 + dVr );

delta_om = theta1 - theta2;


end

%% Local functions

function [ e, h ] = calc_eh( rp, ra, mu )

e = ( ra - rp ) / ( ra + rp ) ;
h = sqrt( 2 .* mu * ra * rp / ( ra + rp) );

end

function [ v_r, v_t ] = vel_LVLH( e, h, theta, mu )

v_t = mu ./ h .* ( 1 + e .* cos(theta) );
v_r = mu ./ h .* e .* sin( theta );

end

