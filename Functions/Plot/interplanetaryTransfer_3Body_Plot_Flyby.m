function [] = interplanetaryTransfer_3Body_Plot_Flyby( planet1, planet2, planet3, departureTime_mjd, flybyTime_mjd, arrivalTime_mjd )

mu_planet2 = astroConstants(10+planet2);
mu_sun = astroConstants(4);

[ r_P1_planet1, ~ ] = uplanet_rv( departureTime_mjd, planet1 ); % Planet1 ephemeris @ P1
[ r_P2_planet2, V_planet2 ] = uplanet_rv( flybyTime_mjd, planet2 );     % Planet2 ephemeris @ P2
[ r_P3_planet3, ~ ] = uplanet_rv( arrivalTime_mjd, planet3 );   % Planet3 ephemeris @ P3

TOF1 = ( flybyTime_mjd - departureTime_mjd ) * 24 * 60 * 60;               % [s]
TOF2 = ( arrivalTime_mjd - flybyTime_mjd) * 24 * 60 * 60;                  % [s]

% Calculate first Lambert arc
[ ~, ~, ~, ~, ~, V_sc_f_1, ~, ~ ] = lambertMR( r_P1_planet1, r_P2_planet2, TOF1, mu_sun, 0, 0, 0, 1 );
% Calculate second Lambert arc
[ ~, ~, ~, ~, V_sc_i_2, ~, ~, ~ ] = lambertMR( r_P2_planet2, r_P3_planet3, TOF2, mu_sun, 0, 0, 0, 1 );

v_inf_minus = V_sc_f_1' - V_planet2; % All vectors in HECI frame
v_inf_m = norm(v_inf_minus);
v_inf_plus = V_sc_i_2' - V_planet2; % All vectors in HECI frame
v_inf_p = norm(v_inf_plus);
u = cross(v_inf_minus,v_inf_plus) ./ cross(v_inf_minus,v_inf_plus);% Direction around which v_inf rotates

[dV_poweredGA,rp,delta,beta_minus,d_m] = poweredGA(  planet2, v_inf_minus, v_inf_plus );

% Find rp vector

rp_vect = rotate_vector( rp * v_inf_minus./v_inf_m, u, -beta_minus );

% Find deltaV
v_P_minus = sqrt( v_inf_m.^2 + 2 .* mu_planet2 ./ rp );
v_P_plus = sqrt( v_inf_p.^2 + 2 .* mu_planet2 ./ rp );
dV_poweredGA = abs( v_P_plus - v_P_minus); % Given by satellite

% Find velocity vectors at pericenter
v_P_unit = rotate_vector(v_inf_minus./v_inf_m, u , d_m/2);
v_P_minus_vect = v_P_minus * v_P_unit;
v_P_plus_vect = v_P_plus * v_P_unit;

t_end = 7000;
points = 1000;
TSPAN_m = linspace(0,-t_end,points);
TSPAN_p = linspace(0,t_end,points);
perturbationModel.type = 'unperturbed';
relTol = 1e-16;
absTol = 1e-16;

[ r_m, ~, ~ ] = orbit_propagation( rp_vect, v_P_minus_vect, TSPAN_m, mu_planet2, relTol, absTol, perturbationModel );
[ r_p, ~, ~ ] = orbit_propagation( rp_vect, v_P_plus_vect, TSPAN_p, mu_planet2, relTol, absTol, perturbationModel );

% Plot results
planet3D( planet2, [0,0,0], 1 );
hold on
grid on

linestyle = '-';
figureNum = gcf;
linewidth = 2.5;
LegendVisibility = 1;

orbitPlotrv( r_m , linestyle, figureNum, [0,0,1], linewidth, LegendVisibility);
orbitPlotrv( r_p , linestyle, figureNum, [1,0,0], linewidth, LegendVisibility);

Delta_m = rp * sqrt( 1 + 2 * mu_planet2 ./ ( rp * norm(v_inf_minus).^2) );
Delta_p = rp * sqrt( 1 + 2 * mu_planet2 ./ ( rp * norm(v_inf_plus).^2) );

Delta_m_vect = rotate_vector( Delta_m * v_inf_minus ./ norm(v_inf_minus), u, -pi/2);
Delta_p_vect = rotate_vector( Delta_p * v_inf_plus ./ norm(v_inf_plus), u, -pi/2);

asymp_m = v_inf_minus./norm(v_inf_minus) * [-50000, 16000] + Delta_m_vect;
asymp_p = v_inf_plus./norm(v_inf_plus) * [-5000, 70000] + Delta_p_vect;

plot3(asymp_m(1,:), asymp_m(2,:), asymp_m(3,:), 'linewidth', 2.5)
plot3(asymp_p(1,:), asymp_p(2,:), asymp_p(3,:), 'linewidth', 2.5)

legend('Incoming leg','Outcoming leg')

hold off

end

%% Local functions

% Powered gravity assist @ planet2
function [dV_poweredGA,rp,delta,beta_minus,d_m] = poweredGA(  planet2, v_inf_minus_vect, v_inf_plus_vect )

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

    beta_minus = deg2rad(90) - d_m(rp)/2;
    d_m = d_m(rp);
end

end
