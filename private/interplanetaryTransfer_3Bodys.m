function [ dV_sat, dV_planet1, dV_planet3, dV_poweredGA, TOF1, TOF2 ] = interplanetaryTransfer_3Bodys( departureTime_mjd, flybyTime_mjd, arrivalTime_mjd, planet1, planet2, planet3 )
% Function calculates total dV required for interplanetary transfer from
% planet1 to planet3, with an intermediate powered gravity-assist at
% planet2.
% Planetocentric phases are neglected -> S/C at initial and final planets is
% assumed to have the same velocity as the planets themselves (on the
% heliocentric orbits). This means that dV_planet = norm( V_planet - V_lambertarc )
% 
% Phase 1: Lambert arc from planet 1 to planet 2 (heliocentric)
%          Goes from departureDAte to flybyDate
% Phase 2: Powered gravity-assist (fly-by) around planet 2 (planet-centered)
%          Seen as impulsive maneuver from heliocentric POV
% Phase 3: Lambert arc from planet 2 to planet 3 (heliocentric)
%          Goes from flybyDate to arrivalDate
%
% P1 conditions @ departureDate -> P2 conditions @ flybyDate  -> P3 conditions @ arrivalDate
%
% ----------------------------------------------------------------------------------------------
%
% PROTOTYPE:
%   [ dV_sat, dV_planet1, dV_planet3, dV_poweredGA, TOF1, TOF2 ] = interplanetaryTransfer_3Body( departureTime_mjd, flybyTime_mjd, arrivalTime_mjd, planet1, planet2, planet3 )
%
% INPUTS:
%   departureTime_mjd[1]  [days]     mjd2000 time at departure
%   flybyTime_mjd[1]      [days]     mjd2000 time at flyby
%   arrivalTime_mjd[1]    [days]     mjd2000 time at arrival
%   planet1[1]            [-]        Integer identifying the first celestial body around the sun
%   planet2[1]            [-]        Integer identifying the second celestial body around the sun
%   planet3[1]            [-]        Integer identifying the third celestial body around the sun
%                                          1:   Mercury
%                                          2:   Venus
%                                          3:   Earth
%                                          4:   Mars
%                                          5:   Jupiter
%                                          6:   Saturn
%                                          7:   Uranus
%                                          8:   Neptune
%                                          9:   Pluto
%                                          10:  Sun
%
% OUTPUTS:
%   dV_sat[1]             [km/s]           Total dV given ONLY by s/c (no natural dV) -> dV_sat = dV_GA + dV_planet1 + dV_planet3
%   dV_planet1[1]         [km/s]           dV given by s/c at planet1 to move from initial heliocentric (planet1) orbit onto first Lambert arc
%   dV_planet3[1]         [km/s]           dV given by s/c at planet3 to move from second Lambert arc onto final heliocentric (planet3) orbit
%   dV_poweredGA[1]       [km/s]           dV given ONLY by s/c at pericenter of planet2 fly-by orbit 
%   TOF1[1]               [s]              TOF on first Lambert arc
%   TOF2[1]               [s]              TOF on second Lambert arc
%
% ----------------------------------------------------------------------------------------------

% Planetary constants
mu_sun = astroConstants(4);

% Calculate planet ephemeris at P1, P2, P3
[ r_P1_planet1, V_P1_planet1 ] = uplanet_rv( departureTime_mjd, planet1 ); % Planet1 ephemeris @ P1
[ r_P2_planet2, V_P2_planet2 ] = uplanet_rv( flybyTime_mjd, planet2 );     % Planet2 ephemeris @ P2
[ r_P3_planet3, V_P3_planet3 ] = uplanet_rv( arrivalTime_mjd, planet3 );   % Planet3 ephemeris @ P3

% Calculate TOF on first and second Lambert arc
TOF1 = ( flybyTime_mjd - departureTime_mjd ) * 24 * 60 * 60;               % [s]
TOF2 = ( arrivalTime_mjd - flybyTime_mjd) * 24 * 60 * 60;                  % [s]

% Set limitations on TOF:
% TOF cannot be negative
if TOF1 < 0 || TOF2 < 0
    dV_sat = NaN; dV_planet1 = NaN; dV_planet3 = NaN; dV_poweredGA = NaN; TOF1 = NaN; TOF2 = NaN;
    return;

end

% % lambertMR.m parameters
% orbit_type = 0;   % Counter-clockwise (prograde) transfer
% Nrev = 0;         % Number of revolutions
% Ncase = 0;        % Small-a solution (not necessary for Nrev = 0)
% optionsLMR = 1;   % Display warnings for convergence

% Calculate first Lambert arc
[ ~, ~, ~, ~, V_sc_i_1, V_sc_f_1, ~, ~ ] = lambertMR( r_P1_planet1, r_P2_planet2, TOF1, mu_sun, 0, 0, 0, 1 );
% Calculate second Lambert arc
[ ~, ~, ~, ~, V_sc_i_2, V_sc_f_2, ~, ~ ] = lambertMR( r_P2_planet2, r_P3_planet3, TOF2, mu_sun, 0, 0, 0, 1 );

% Calculate dV given at powered GA with poweredGA()
v_inf_minus_vect_P2 = V_sc_f_1' - V_P2_planet2;
v_inf_plus_vect_P2 = V_sc_i_2' - V_P2_planet2;
dV_poweredGA = poweredGA(  planet2, v_inf_minus_vect_P2, v_inf_plus_vect_P2 );

if isnan(dV_poweredGA)

    % Maneuver is below minimum acceptable rp at flyby planet
    dV_sat = NaN; dV_planet1 = NaN; dV_planet3 = NaN; dV_poweredGA = NaN; TOF1 = NaN; TOF2 = NaN;

else

    % Calculate dV given at planet1

    dV_planet1 = norm( V_sc_i_1' - V_P1_planet1 );

%     v_inf_plus_vect_P1 = V_sc_i_1' - V_P1_planet1;
%     Delta_hyp_planet1 = r_park1 .* sqrt( 1 + 2 .* mu_planet1 ./ ( r_park1 .* norm( v_inf_plus_vect_P1 ).^2 ) );
%     h_hyp_planet1 = Delta_hyp_planet1 .* norm( v_inf_plus_vect_P1 );
%     v_P_hyp_planet1 = h_hyp_planet1 ./ r_park1; % Velocity at pericenter of hyperbola @ planet1
%     v_c_park1 = sqrt( mu_planet1 ./ r_park1); % Velocity at pericenter of circular orbit @ planet1
%     dV_park1 = abs( v_P_hyp_planet1 - v_c_park1 );

    % Calculate dV given at planet3
    dV_planet3 = norm( V_sc_f_2' - V_P3_planet3 );



%     v_inf_minus_vect_P3 = V_sc_f_2' - V_P3_planet3;
%     Delta_hyp_planet3 = r_park3 .* sqrt( 1 + 2 .* mu_planet3 ./ ( r_park3 .* norm( v_inf_minus_vect_P3 ).^2 ) );
%     h_hyp_planet3 = Delta_hyp_planet3 .* norm( v_inf_minus_vect_P3 );
%     v_P_hyp_planet3 = h_hyp_planet3 ./ r_park3; % Velocity at pericenter of hyperbola @ planet3
%     v_c_park3 = sqrt( mu_planet3 ./ r_park3); % Velocity at pericenter of circular orbit @ planet3
%     dV_park3 = abs( v_P_hyp_planet3 - v_c_park3 );

    % dV given ONLY by satellite ( excluding natural dV given by gravity-assist )
    dV_sat = dV_poweredGA + dV_planet1 + dV_planet3;

    % TOTAL dV ( includes dV given by satellite and dV given naturally by gravity-assist )
    %dV_tot = norm( V_sc_f_1 - V_sc_i_2 ) + dV_planet1 + dV_planet3;

end

    if isnan(dV_sat)
        dV_sat = 1e10;
    end

end



%% Local functions

% Powered gravity assist @ planet2
function dV_poweredGA = poweredGA(  planet2, v_inf_minus_vect, v_inf_plus_vect )

mu_planet2 = astroConstants( 10+planet2 );
R_planet2 = astroConstants( 20+planet2 );

% Norms of v_inf vectors
v_inf_m = norm(v_inf_minus_vect);
v_inf_p = norm(v_inf_plus_vect);

% Calculate turning angle for powered gravity-assist
delta = acos( dot( v_inf_minus_vect, v_inf_plus_vect ) ./ ( v_inf_m .* v_inf_p ) );

% Calculate non-linear function FUN = FUN(rp) = 0 to solve for rp numerically
e_m = @(rp) 1 + rp .* v_inf_m.^2 ./ mu_planet2;
d_m = @(rp) 2 .* asin( 1 ./ e_m(rp) );
e_p = @(rp) 1 + rp .* v_inf_p.^2 ./ mu_planet2;
d_p = @(rp) 2 .* asin( 1 ./ e_p(rp) );

FUN = @(rp) delta - d_m(rp) ./ 2 - d_p(rp) ./ 2;

% Set lowest acceptable height of pericenter at the flyby
h_atm = 100;                        % Altitude of atmosphere of planet2
rp_min = R_planet2 + h_atm;         % Minimum acceptable rp for powered GA

% Find zero of FUN corresponding to rp
rp_GUESS = rp_min;

if ~isinf(FUN(rp_GUESS)) & isreal(FUN(rp_GUESS))

    options = optimset('Display','off');
    rp = fzero( FUN, rp_GUESS, options );

else

    dV_sat = NaN; dV_planet1 = NaN; dV_planet3 = NaN; dV_poweredGA = NaN; TOF1 = NaN; TOF2 = NaN;
    return;

end


%Check validity of calculated rp ( cannot be smaller than rp_minimum )
if rp < rp_min

    % Condition is not verified, set dV_poweredGA = NaN and exit
    dV_poweredGA = NaN;

else
    
    % Impact parameters of incoming and outgoing hyperbolas
    Delta_m = rp .* sqrt(1 + 2 .* mu_planet2 ./ ( rp .* v_inf_m.^2 ) );
    Delta_p = rp .* sqrt(1 + 2 .* mu_planet2 ./ ( rp .* v_inf_p.^2 ) );

    % Angular momentums of incoming and outgoing hyperbolas
    h_hyp_m = Delta_m .* v_inf_m;
    h_hyp_p = Delta_p .* v_inf_p;

    % Calculate dV at pericenter ( given by satellite --> does not include natural dV given by planet2)
    v_P_m = h_hyp_m ./ rp;
    v_P_p = h_hyp_p ./ rp;

    dV_poweredGA = abs( v_P_m - v_P_p );

    % Planet of hyperbola is defined if i know v_inf- and v_inf+
    % u = cross(v_inf-,v_inf+) ./ norm( cross(v_inf-,v_inf+) ) -> Normal to
    % plane

end

end