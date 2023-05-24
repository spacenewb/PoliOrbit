function [ dV_tot, a2, e2, t_transfer ] = phasing_maneuver_CIRCULAR( dTheta, a1,0, i1, OM1, om1, mu, n_rev, R_planet )
% ONLY WORKS FOR A CIRCULAR ORBIT
% phasing_maneuver calculates necessary dV and orbit to obtain an imposed 
% dTheta = theta_final - theta_initial. This function considers as
% maneuvering point the pericenter of the orbit, and a number of
% revolutions on transfer orbit before complete phasing n_rev.
%
% PROTOTYPE:
%   [ dV_tot, a2, e2, t_transfer ] = phasing_maneuver( dTheta, a1, e1, i1, OM1, om1, mu, n_rev, R_planet )
%
% INPUTS:
%   dTheta    [1]    Final-Initial real anomaly on orbit [rad]
%   a1        [1]    Semi-major axis                     [km]
%   e1        [1]    Eccentricity                        [-]
%   i1        [1]    Inclination                         [rad]
%   OM1       [1]    RAAN                                [rad]
%   om1       [1]    Argument of pericenter              [rad]
%   mu        [1]    Planetary constant                  [km^3/s^2]
%   n_rev     [1]    Number of revolutions               [-]
%   R_planet  [1]    If specified, mean radius of planet [km]
%
% OUTPUTS:
%   dV_tot    [1]    Total deltaV required               [km/s]
%   a_t       [1]    Semi-major axis of transfer orbit   [km]
%   e_t       [1]    Eccentricity of transfer orbit      [-]
%   t_t       [1]    Time to complete phasing            [s]  
% Other parameters of orbit remain unchanged
%   
% CONTRIBUTORS:
%   Matteo D'Ambrosio

[ T1, rp1, ra1, ~, ~ ] = orbitParameters( a1, e1, i1, OM1, om1, 0, mu );
h1 = sqrt( 2 .* mu .* ra1 .* rp1 ./ ( ra1 + rp1) );

% Transfer orbit
T_T = T1 + T1 .* ( -sign(dTheta) ) .* abs(dTheta) ./ ( 2 .* pi .* n_rev ); % Period of tranfer orbit: I have to increase or decrease the period on the orbit based on the dTheta i need.
a2 = ( mu .* T_T.^2 ./ ( 4 .* pi.^2 ) ).^( 1/3 );
t_transfer = n_rev * T_T;

switch sign(dTheta)

    case 1 % T_T < T1
        
        ra2 = rp1;
        rp2 = 2 .* a2 - ra2;
        e2 = abs((ra2 - rp2) ./ (ra2 + rp2));
        h2 = sqrt( 2 .* mu .* ra2 .* rp2 ./ ( ra2 + rp2) );

        dV_tot = 2 .* abs( ( h2 ./ ra2 - h1 ./ rp1) );

    case -1 % T_T > T1

        rp2 = rp1;
        ra2 = 2 .* a2 - rp2;
        e2 = abs((ra2 - rp2) ./ (ra2 + rp2));
        h2 = sqrt( 2 .* mu .* ra2 .* rp2 ./ ( ra2 + rp2) );

        dV_tot = 2 .* abs( ( h2 ./ rp2 - h1 ./ rp1 ) );

end

if nargin == 9
    if rp2 < R_planet || ra2 < R_planet

        error('This maneuver impacts the planet. Consider dividing dTheta by 2 and waiting for 2 revolutions on the orbit')

    end
end

end